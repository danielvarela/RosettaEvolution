
#include <map>
#include <vector>
#include <string>
#include <stdlib.h>     /* exit, EXIT_FAILURE **/

#include "DEoperator.hpp"

DE_Operator::DE_Operator() {
  init_protein_archive();

  // score per stage
  score_per_stage["stage1"] = "score0";
  score_per_stage["stage2"] = "score1";
  score_per_stage["stage3"] = "score2";
  score_per_stage["stage4"] = "score3";

  // gmax_per_stage;
  gmax_per_stage["stage1"] = 30;
  gmax_per_stage["stage2"] = 10;
  gmax_per_stage["stage3"] = 190;
  gmax_per_stage["stage4"] = 100000;

  // gmax_per_stage["stage1"] = 30;
  // gmax_per_stage["stage2"] = 10;
  // gmax_per_stage["stage3"] = 190;
  // gmax_per_stage["stage4"] = 100000;

  distances_map["rmsd"] = rmsd;
  distances_map["rmsd_native_diff"] = rmsd_native_diff;
  distances_map["euclidean"] = euclidean;
  distances_map["euclidean_loop"] = euclidean_loop;
  distances_map["euclidean_diff_abs"] = euclidean_diff_abs;
  distances_map["euclidean_mario"] = euclidean_mario;
  distances_map["euclidean_partial_mario"] = euclidean_partial_mario;

  protocol_name_map["Shared"] = Shared;
  protocol_name_map["HybridShared"] = HybridShared;
  protocol_name_map["CrowdingDE"] = CrowdingDE;
  protocol_name_map["HybridCrowdingDE"] = HybridCrowdingDE;

  init_popul_strategy_map["total_random"] = total_random;
  init_popul_strategy_map["random_pose_based"] = random_pose_based;
  init_popul_strategy_map["total_random_pose_based"] = total_random_pose_based;
  init_popul_strategy_map["init_popul_with_stage"] = init_popul_with_stage;

  fragment_insertion_strategy_map["my_frag_insertion"] = my_frag_insertion;
  fragment_insertion_strategy_map["greedy_search"] = greedy_search;
  fragment_insertion_strategy_map["stage_rosetta_mover"] = stage_rosetta_mover;
  fragment_insertion_strategy_map["ILS_as_julia"] = ILS_as_julia;

}


void
DE_Operator::init_available_stages() {
  boost::char_separator<char> sep(",");
  std::string input_stages = app_options.get<std::string>("Protocol.stages");
  trim(input_stages);
  boost::tokenizer<boost::char_separator<char>> tokens(input_stages, sep);
  for (boost::tokenizer<boost::char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it) {
    vec_stages.push_back(*it);
  }
}

DE_Operator::DE_Operator(std::string prot_selected, boost::property_tree::ptree opt) : DE_Operator() {
  app_options = opt;
  fit_radius = app_options.get<double>("Extra.fitrad");
  init_files(prot_selection[prot_selected]);
  init_setup();
  init_two_stages_mover();
  init_available_stages();
}

DE_Operator::DE_Operator(std::string prot_selected) : DE_Operator() {
  init_files(prot_selection[prot_selected]);
  init_setup();
  init_default_score();
  init_popul = initialize_init_popul_strategy(std::string("total_random_pose_based"));
}


void
DE_Operator::init_protein_archive() {
  Protinfo sap;
  sap.name = "1sap";
  sap.pdb_file = "./1sap_A.pdb";
  sap.map_file = "./4054.map";
  sap.frag_3 = "./frag3";
  sap.frag_9 = "./frag9";
  sap.ss_file = "./1sapA.ss2";

  Protinfo j9s;
  j9s.name = "3j9s";
  j9s.pdb_file = "./input.pdb";
  j9s.map_file = "./emd_6272.map";
  j9s.ss_file = "./sec.ss2";

  Protinfo c8c;
  c8c.name = "1c8cA";
  c8c.pdb_file = "./input_files/info_1c8cA/vf_1c8c.pdb";
  c8c.map_file = "./emd_6272.map";
  c8c.frag_3 = "./input_files/info_1c8cA/boinc_vf_aa1c8cA03_05.200_v1_3";
  c8c.frag_9 = "./input_files/info_1c8cA/boinc_vf_aa1c8cA09_05.200_v1_3";
  c8c.ss_file = "./input_files/info_1c8cA/vf_1c8cA.psipred_ss2";

  Protinfo elw;
  elw.name = "1elwA";
  elw.pdb_file = "./input_files/info_1elwA/vf_1elw.pdb";
  elw.map_file = "./emd_6272.map";
  elw.frag_3 = "./input_files/info_1elwA/boinc_vf_aa1elwA03_05.200_v1_3";
  elw.frag_9 = "./input_files/info_1elwA/boinc_vf_aa1elwA09_05.200_v1_3";
  elw.ss_file = "./input_files/info_1elwA/vf_1elwA.psipred_ss2";

  Protinfo opd;
  opd.name = "1opd";
  opd.pdb_file = "./input_files/info_1opd/vf_1opd.pdb";
  opd.map_file = "./emd_6272.map";
  opd.frag_3 = "./input_files/info_1opd/boinc_vf_aa1opd_03_05.200_v1_3";
  opd.frag_9 = "./input_files/info_1opd/boinc_vf_aa1opd_09_05.200_v1_3";
  opd.ss_file = "./input_files/info_1opd/vf_1opd_.psipred_ss2";

  prot_selection["1c8cA"] = c8c;
  prot_selection["1elwA"] = elw;
  prot_selection["1opd"] = opd;
}

void
DE_Operator::init_setup() {
  pose_size = ipose.size();
  std::cout << "pose size " << pose_size << std::endl;
  pose_ = core::pose::PoseOP(new core::pose::Pose(ipose));
  best_pose_ = pose_->clone();
  native_pose_ = pose_->clone();
  scorefxn = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score3").c_str());
  std::cout << "start score: " << (*scorefxn)(*pose_) << std::endl;
  de_len = 10;
  start_res = 2;
  best_score = 1000000;
  // COMPLETE ABINITIO MOVER AT THE BOTTOM OF THE FILE
}

CalculateDistancePopulationPtr
DE_Operator::use_distances_strategy(std::string option) {
  switch (distances_map[option]) {
  case rmsd: {
    calculate_distances_popul = CalculateDistancePopulationPtr( new CalculateRmsdDistancePopulation(native_pose_, ffxn, ss, scorefxn, fit_radius));
    break;
  }
  case rmsd_native_diff: {
    calculate_distances_popul = CalculateDistancePopulationPtr( new CalculateNativeDiffDistancePopulation(native_pose_, ffxn, ss, scorefxn, fit_radius));
    break;
  }
  case euclidean: {
    calculate_distances_popul = CalculateDistancePopulationPtr( new CalculateEuclideanDistancePopulation(native_pose_, ffxn, ss, scorefxn, fit_radius));
    break;
  }
  case euclidean_loop: {
    calculate_distances_popul = CalculateDistancePopulationPtr( new CalculateEuclideanLoopDistancePopulation(native_pose_, ffxn, ss, scorefxn, fit_radius));
    break;
  }
  case euclidean_diff_abs: {
    calculate_distances_popul = CalculateDistancePopulationPtr( new CalculateEuclideanDiffAbsDistancePopulation(native_pose_, ffxn, ss, scorefxn, fit_radius));
    break;
  }
  case euclidean_partial_mario: {
    calculate_distances_popul = CalculateDistancePopulationPtr( new CalculateEuclideanMarioPartialDistancePopulation(native_pose_, ffxn, ss, scorefxn, fit_radius));
    break;
  }
  case euclidean_mario: {
    calculate_distances_popul = CalculateDistancePopulationPtr( new CalculateEuclideanMarioDistancePopulation(native_pose_, ffxn, ss, scorefxn, fit_radius));
    break;
  }
  default:
    std::cout << "error with distances strategy" << std::endl;
    exit(1);
    break;
  }
  return calculate_distances_popul;
}


void
DE_Operator::init_files(Protinfo prot_info) {
  std::cout << "init protein " << prot_info.name << std::endl;
  read_pose(prot_info.pdb_file, ipose);
  ss = read_ss2(prot_info.ss_file, ipose);
  std::cout << "start read fragment lib " << std::endl;
  read_frag_lib(prot_info.frag_3, frag_set_ );
  read_frag_lib(prot_info.frag_9, frag_set_large );
  frag_opt.frag_set_ = frag_set_;
  frag_opt.frag_set_large = frag_set_large;
  //    read_density_map(prot_info.map_file, ipose);
}

boost::shared_ptr<MoverDE>
DE_Operator::prepare_stage(std::string stage_name) {
  /* INIT all funcions and classes neeeded by Differential Evolution
   */
  scorefxn = core::scoring::ScoreFunctionFactory::create_score_function(std::string(score_per_stage[stage_name]).c_str());
  frag_opt.scorefxn = scorefxn;
  frag_opt.stage_name = stage_name;

  bool fragment_at_trials = false;
  boost::property_tree::ptree::const_assoc_iterator it = app_options.find("Fragments.strategy_at_trials");
  if(app_options.not_found() == it) {
    fragment_at_trials = true;
    frag_mover = initialize_fragment_insertion_strategy(app_options.get<std::string>("Fragments.strategy_at_trials"));
  } else {
    fragment_at_trials = false;
    frag_mover = initialize_fragment_insertion_strategy("my_frag_insertion");
  }

   if (fragment_at_trials) {
     ffxn = boost::shared_ptr<FitFunction>( new PoseFragmentFunction(pose_, scorefxn, ss, frag_mover));

  } else {
    ffxn = boost::shared_ptr<FitFunction>( new PoseScoreFunction(pose_, scorefxn, ss, frag_mover));
  }


   if ( (app_options.get<std::string>("Protocol.name") =="HybridShared") ||
	(app_options.get<std::string>("Protocol.name") =="HybridCrowdingDE")) {
     bool fragment_at_popul = false;
     boost::property_tree::ptree::const_assoc_iterator it = app_options.find("Fragments.strategy_at_population");
     if(app_options.not_found() == it) {
       fragment_at_popul = true;
       initialize_local_search_to_apply_at_population(app_options.get<std::string>("Fragments.strategy_at_population"));
     } else {
       //default option
       initialize_local_search_to_apply_at_population(std::string("stage_rosetta_mover"));
     }
   }


  ffxn->name_ = std::string("Score Function for stage " + stage_name);
  print_best = PrintBestIndividualPtr( new PrintBestIndividual(pose_, ffxn, ss, scorefxn));

  //calculate_distances_popul = use_distances_strategy("rmsd");
  calculate_distances_popul = use_distances_strategy( app_options.get<std::string>("Protocol.distance_strategy") );

  update_current_population_to_stage_score();

  boost::shared_ptr<MoverDE> de = init_differential_evolution_protocol();
  de->Gmax = gmax_per_stage[stage_name];
  de->use_print_class = true;
  de->print_best_ind = print_best;
  de->calculate_distances_popul = calculate_distances_popul;

  //  de->popul = current_population;

  return de;
}

void
DE_Operator::initialize_local_search_to_apply_at_population(std::string input_option) {
  local_search_frag_mover = initialize_fragment_insertion_strategy(input_option);
  local_search = boost::shared_ptr<LocalSearchIndividualMover>(new LocalSearchIndividualMover(pose_, scorefxn, ss, local_search_frag_mover));
}

boost::shared_ptr<FragInsertionMover>
DE_Operator::initialize_fragment_insertion_strategy(std::string input_option ) {
 boost::shared_ptr<FragInsertionMover> strategy_return;

  switch (fragment_insertion_strategy_map[input_option]) {
  case my_frag_insertion: {
    strategy_return = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::my_frag_insertion, frag_opt);
    break;
  }
  case greedy_search: {
    strategy_return = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::greedy_search, frag_opt);
    break;
  }
  case stage_rosetta_mover: {
    strategy_return = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::stage_rosetta_mover, frag_opt);
    break;
  }
  case ILS_as_julia: {
    strategy_return = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::ILS_as_julia, frag_opt);
    break;
  }
  default:
    std::cout << "no insertion strategy option" << std::endl;
    exit(1);
    break;
  }
  return strategy_return;


}

boost::shared_ptr<InitPopulation>
DE_Operator::initialize_init_popul_strategy(std::string init_popul_option) {
  boost::shared_ptr<InitPopulation> init_popul_in_options_file;

  switch (init_popul_strategy_map[init_popul_option]) {
  case total_random: {
    init_popul_in_options_file = boost::shared_ptr<InitPopulation>(new InitPopulation(ffxn));
    break;
  }
  case total_random_pose_based: {
    init_popul_in_options_file = boost::shared_ptr<InitPopulation>(new PoseTotalRandomInitPopulation(ffxn, pose_, ss));
    break;
  }
  case random_pose_based: {
    init_popul_in_options_file = boost::shared_ptr<InitPopulation>(new PosePopulation(ffxn, pose_, ss));
    break;
  }
  case init_popul_with_stage: {
    init_popul_in_options_file = boost::shared_ptr<InitPopulation>(new TwoStagesInitPopulation(ffxn, pose_, two_stages_mover, ss));
    break;
  }
  default:
    std::cout << "problem selecting the init population strategy" << std::endl;
    exit(1);
    break;
  }
  return init_popul_in_options_file;
}

void
DE_Operator::init_default_score() {
  core::scoring::ScoreFunctionOP score3 = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score3").c_str());
  frag_opt.scorefxn = score3;
  frag_opt.stage_name = "stage4";
  frag_mover = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::greedy_search, frag_opt);
  ffxn = boost::shared_ptr<FitFunction>( new PoseScoreFunction(pose_, score3, ss, frag_mover)); 
}

void
DE_Operator::init_two_stages_mover() {
  core::scoring::ScoreFunctionOP score1 = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score1").c_str());
  frag_opt.scorefxn = score1;
  frag_opt.stage_name = "stage1";
  frag_mover = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::greedy_search, frag_opt);
  ffxn = boost::shared_ptr<FitFunction>( new PoseScoreFunction(pose_, score1, ss, frag_mover));
  two_stages_mover = boost::shared_ptr<InitStagesMover>( new InitStagesMover(score1, frag_set_, frag_set_large));
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
    for (int i = 0; i < current_population.size(); i++) {
      ffxn->score(current_population[i]);
    }
  }
}

boost::shared_ptr<MoverDE>
DE_Operator::init_differential_evolution_protocol() {
  //     boost::shared_ptr<MoverDE> de = boost::shared_ptr<MoverDE>(new SharedMoverDE(ffxn, init_popul));
  ConfigurationDE conf = ConfigurationDE(app_options.get<int>("DE.NP"), app_options.get<double>("DE.CR"), app_options.get<double>("DE.F"), app_options.get<double>("Extra.fitrad") , app_options.get<std::string>("Protocol.prot"));

  std::string protocol_name = app_options.get<std::string>("Protocol.name");

  boost::shared_ptr<MoverDE> de ;
  switch (protocol_name_map[protocol_name]) {
  case Shared: {
    de = boost::shared_ptr<MoverDE>(new SharedMoverDE(app_options, ffxn, current_population));
    break;
  }
  case HybridShared: {
    de = boost::shared_ptr<MoverDE>(new SharedHybridMoverDE(app_options, ffxn, current_population, local_search));
    break;
  }
  case CrowdingDE: {
    de = boost::shared_ptr<MoverDE>(new CrowdingMoverDE(app_options, ffxn, current_population));
    break;
  }
  case HybridCrowdingDE: {
    de = boost::shared_ptr<MoverDE>(new CrowdingHybridMoverDE(app_options, ffxn, current_population, local_search));
    break;
  }
default:
  std::cout << "error with Protocol Name" << std::endl;
  exit(1);
    break;
  }
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
DE_Operator::run() {
  std::string stage_name;
  std::cout << "run " << std::endl;

  initialize_population_stage1();

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
  }

  std::cout << "================================================" << std::endl;
  print_best_pose();
}
