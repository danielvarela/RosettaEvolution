
#include <map>
#include <vector>
#include <string>

#include <stdlib.h>     /* exit, EXIT_FAILURE **/

#include "DE_Operator.hpp"
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>

#include "../MpiFiles/MasterRosettaCalculator.hpp"
#include "ProteinInfoArchive.hpp"
#include "FragmentInserterBuilder.hpp"
#include "StageBuilder.hpp"
#include "OptionsMapInitializer.hpp"
#include "PopulStrategyInitializer.hpp"
#include "FitnessFunctionInitializer.hpp"

DE_Operator::DE_Operator() {
  init_protein_archive();
  score_per_stage = OptionsMapInitializer::score_per_stage_build();
  gmax_per_stage = OptionsMapInitializer::gmax_per_stage_build();
  distances_map = OptionsMapInitializer::distances_map_build();
  protocol_name_map = OptionsMapInitializer::protocol_name_map_build();
  init_popul_strategy_map = OptionsMapInitializer::init_popul_strategy_map_build();
}

void
DE_Operator::init_available_stages() {
  boost::char_separator<char> sep(",");
  std::string input_stages = app_options.get<std::string>("Protocol.stages");
  vec_stages = OptionsMapInitializer::parse_stages(input_stages);
}

DE_Operator::DE_Operator(std::string prot_selected, boost::property_tree::ptree opt) : DE_Operator() {
  app_options = opt;
  fit_radius = app_options.get<double>("Extra.fitrad");
  init_files(prot_selection[prot_selected]);
  init_setup();
  init_defaults();
  init_two_stages_mover();
  init_available_stages();
}

DE_Operator::DE_Operator(std::string prot_selected) : DE_Operator() {
  init_files(prot_selection[prot_selected]);
  init_setup();
  init_defaults();
  init_popul = initialize_init_popul_strategy(std::string("total_random_pose_based"));
}

void
DE_Operator::init_protein_archive() {
  ProteinInfoArchive prot_archive;
  prot_selection = prot_archive.get_map();
}

void
DE_Operator::init_setup() {
  pose_size = ipose.size();
  pose_ = core::pose::PoseOP(new core::pose::Pose(ipose));
  best_pose_ = pose_->clone();
  native_pose_ = pose_->clone();
  scorefxn = OptionsMapInitializer::CreateScoreFunctionStrategy::get()->build(std::string("score3"));
  //StageBuilder stage_builder(app_options, scorefxn, frag_opt);
  //stage_builder.initializer_distances_calculator();
  // calculate_distances_popul = stage_builder.calculate_distances_popul;
  calculate_distances_popul =use_distances_strategy("rmsd");
  de_len = 10;
  start_res = 2;
  best_score = 1000000;
}

CalculateDistancePopulationPtr
DE_Operator::use_distances_strategy(std::string option) {
  return CalculateRmsdDistancePopulation(native_pose_, ffxn, ss, scorefxn, fit_radius).use_distances_strategy(option);
}

void
DE_Operator::init_files(Protinfo prot_info) {
  //std::cout << "init protein " << prot_info.name << std::endl;
  read_pose(prot_info.pdb_file, ipose);
  ss = read_ss2(prot_info.ss_file, ipose);
  //std::cout << "start read fragment lib " << std::endl;
  frag_opt = frag_inserter_builder.init_frag_files(prot_info, app_options);

}

boost::shared_ptr<MoverDE>
DE_Operator::prepare_stage(std::string stage_name) {
  /* INIT all funcions and classes neeeded by Differential Evolution
   */
  StageBuilder stage_builder(app_options, scorefxn, frag_opt);
  calculate_distances_popul = use_distances_strategy("rmsd");
  boost::shared_ptr<MoverDE> de = stage_builder.prepare_DE_for_stage( stage_name, current_population);
  return de;
}

void
DE_Operator::initialize_local_search_to_apply_at_population(std::string input_option) {
  local_search_frag_mover = initialize_fragment_insertion_strategy(input_option);
  local_search = boost::shared_ptr<LocalSearchIndividualMover>(new LocalSearchIndividualMover(pose_, scorefxn, ss, local_search_frag_mover));
}

boost::shared_ptr<FragInsertionMover>
DE_Operator::initialize_fragment_insertion_strategy(std::string input_option ) {
  return frag_inserter_builder.get(input_option, frag_opt);
}

boost::shared_ptr<InitPopulation>
DE_Operator::initialize_init_popul_strategy(std::string init_popul_option) {
  return OptionsMapInitializer::PopulStrategyInitializer(app_options, frag_opt).build();
}

void
DE_Operator::init_defaults() {
  core::scoring::ScoreFunctionOP score3 = OptionsMapInitializer::CreateScoreFunctionStrategy::get()->build(std::string("score3"));
  frag_opt.scorefxn = score3;
  frag_opt.stage_name = "stage4";
  frag_mover = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::hybrid_mover, frag_opt);
  ffxn = OptionsMapInitializer::FitnessFunctionInitializer(frag_opt, frag_mover, false).set_up();
}

void
DE_Operator::init_two_stages_mover() {
  init_popul = initialize_init_popul_strategy(app_options.get<std::string>("Protocol.init_strategy") );
}

void
DE_Operator::update_current_population_to_stage_score( ) {
  if (current_population.size() == 0) {
    //    current_population = de->popul;
    std::cout << "error, the current population was not initialize" << std::endl;
    exit(1);
  } else {
    CalculateDistancePopulation::DistancesResult result =
      calculate_distances_popul->run(current_population);
    std::vector<double> current_rmsd = result.rmsd_to_native;
    // rescore population
    ffxn = OptionsMapInitializer::FitnessFunctionInitializer(frag_opt, frag_mover, false).set_up();
  boost::shared_ptr<MasterRosettaCalculator> mpi_calculator;
#if(MPI_ENABLED)
  mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterRosettaCalculator(ffxn));
#endif


#if(MPI_ENABLED)
  int frags_at_popul = 1;
  current_population = mpi_calculator->run(current_population, frags_at_popul);
#else
  for (int i = 0; i < current_population.size(); i++) {
     ffxn->score(current_population[i]);
   }
#endif


  }
}

boost::shared_ptr<MoverDE>
DE_Operator::init_differential_evolution_protocol(std::string option_protocol = "") {
  boost::shared_ptr<MoverDE> de = AlgorithmBuilder(app_options, ffxn,  current_population,  local_search, frag_opt).build(option_protocol);
  return de;
}

void
DE_Operator::print_best_pose() {
  double obtained_score = (*scorefxn)(*best_pose_);
  std::cout << "obtained score: " << obtained_score << std::endl;
  for (int i = 1; i < static_cast<int>(best_pose_->size()); ++i) {
    std::cout << best_pose_->phi(i) << " " << best_pose_->psi(i) << " ";
  }
  std::cout << std::endl;
  best_pose_->dump_pdb("./result_relax.pdb");
}

void
DE_Operator::apply_differential_evolution_with_stage(std::string stage_name) {
  boost::shared_ptr<MoverDE> de;
  de = prepare_stage(stage_name);
  de->apply();
  current_population = de->popul;
  best_pose_ = pose_->clone();
  ffxn->fill_pose(best_pose_, de->best_ind(), ss);
  std::cout << "BEST POSE FOUND AT STAGE " << stage_name << " " << de->best_ind().score << std::endl;
  std::cout << "finish apply " << stage_name << std::endl;

  // boost::mpi::communicator world;
  // int STOP_TAG = -42;
  // for (int i = 1; i < world.size(); i++) {
  //     world.send(i, 0, STOP_TAG);
  // }
  //exit(0);

}

core::pose::PoseOP
DE_Operator::get_native_pose() {
  return native_pose_;
}

void
DE_Operator::initialize_population_stage1() {
  // current_population is modified using the stage1 of rosetta
  std::cout << "=== init the current population ================ " << std::endl;
  init_popul->apply(current_population, app_options.get<int>("DE.NP"), ffxn->D() );
  std::cout << "=== finish initialization current population === " << std::endl;
}

void
DE_Operator::run_complete_abinitio() {
#if(MPI_ENABLED)
  core::scoring::ScoreFunctionOP score_4 = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score3").c_str());
  boost::shared_ptr<CompleteAbinitioMover> abinitio = boost::shared_ptr<CompleteAbinitioMover>(new CompleteAbinitioMover(score_4, frag_opt.frag_set_, frag_opt.frag_set_large));
  prepare_stage("stage4");
  ffxn = boost::shared_ptr<FitFunction>( new PoseFragmentFunction(pose_, scorefxn, ss, abinitio));
#if(USE_CRYO_EM)
  ffxn = boost::shared_ptr<FitFunction>( new PoseDensityFragmentFunction(pose_, scorefxn, ss, frag_mover));
#endif
  ffxn->set_name("complete");
  boost::shared_ptr<FitFunction> simple_ffxn = boost::shared_ptr<FitFunction>( new PoseScoreFunction(pose_, scorefxn, ss, frag_mover));

#if(USE_CRYO_EM)
  simple_ffxn = boost::shared_ptr<FitFunction>( new PoseDensityFunction(pose_, scorefxn, ss, frag_mover));
#endif

  boost::shared_ptr<MasterRosettaCalculator> mpi_calculator;
  mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterRosettaCalculator(ffxn));
  current_population = mpi_calculator->run(current_population);
  std::cout << "finish rosetta" << std::endl;
  print_final_population(current_population);
#endif
}

void DE_Operator::print_final_population(std::vector<Individual> popul ) {
  calculate_distances_popul =use_distances_strategy("rmsd");
  CalculateDistancePopulation::DistancesResult result =
    calculate_distances_popul->run(popul);
  std::vector<double> current_rmsd = result.rmsd_to_native;

  double acc = 0;
  double best = 100000;
  for (int i = 0; i < app_options.get<int>("DE.NP"); i++) {
    double score =  (-1 * SCORE_ERROR_FIXED) +  popul[i].score;
    acc += score;
    if (score < best) {
      best = score;
    }
  }
  std::cout << "[GEN] 99999 " << best << " " << acc/popul.size() << std::endl;
  std::cout << "[POP] ";
  for (int i = 0; i < app_options.get<int>("DE.NP"); i++) {
    double score = popul[i].score;
    std::cout << (-1 * SCORE_ERROR_FIXED) +  score << " , ";
  }

  std::cout << std::endl;
  std::cout << "[RMSD_NATIVE] ";
  for (int i = 0; i < app_options.get<int>("DE.NP"); i++) {
    double rmsd = current_rmsd[i];
    std::cout << rmsd << " , ";
  }
  std::cout << std::endl;
}

void
DE_Operator::run() {
  std::string stage_name;
  std::cout << "run " << std::endl;
  std::string protocol_name = app_options.get<std::string>("Protocol.name");
  std::string prot_name = app_options.get<std::string>("Protocol.prot");

  if (false) {
    initialize_population_stage1();
    run_complete_abinitio();
  } else {
    initialize_population_stage1();
    PrintRMSDvsMarioAnalisys printer_mario_pairs( prot_name, use_distances_strategy("rmsd"),  use_distances_strategy("euclidean_partial_mario"));
    // std::cout << "================================================" << std::endl;
    // printer_mario_pairs.apply(current_population);
    std::cout << "================================================" << std::endl;

    if ( std::find(vec_stages.begin(), vec_stages.end(), "stage2") != vec_stages.end()) {
      stage_name = "stage2";
      std::cout << "================================================" << std::endl;
      std::cout << "START STAGE " << stage_name << std::endl;
      apply_differential_evolution_with_stage(stage_name);
    }


    if ( std::find(vec_stages.begin(), vec_stages.end(), "stage3") != vec_stages.end()) {
      stage_name = "stage3";
      std::cout << "================================================" << std::endl;
      std::cout << "START STAGE " << stage_name << std::endl;
      apply_differential_evolution_with_stage(stage_name);
    }

    if ( std::find(vec_stages.begin(), vec_stages.end(), "stage4") != vec_stages.end()) {
      stage_name = "stage4";
      std::cout << "================================================" << std::endl;
      std::cout << "START STAGE " << stage_name << std::endl;
      apply_differential_evolution_with_stage(stage_name);
      std::cout << "llega a acabar" << std::endl;
    }

    // std::cout << "================================================" << std::endl;
    // printer_mario_pairs.apply(current_population);
    // std::cout << "================================================" << std::endl;

    // printer_mario_pairs.set_resolution(app_options.get<int>("Protocol.map_res"));
    // std::string tag_dist_degrees = "deg" + app_options.get<std::string>("Protocol.disturb_degrees");
    // printer_mario_pairs.set_extra_tags(tag_dist_degrees);
    // printer_mario_pairs.dump_pdb(current_population);
    //print_best_pose();
    print_final_population(current_population);
  }
}
