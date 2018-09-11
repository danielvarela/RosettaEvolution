
#include "DEoperator.hpp"
#include <map>
#include <vector>
#include <stdlib.h>     /* exit, EXIT_FAILURE **/

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

  distances_map["rmsd"] = rmsd;
  distances_map["euclidean"] = euclidean;
  distances_map["euclidean_loop"] = euclidean_loop;

  protocol_name_map["Shared"] = Shared;
  protocol_name_map["HybridShared"] = HybridShared;
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
  init_available_stages();
}

DE_Operator::DE_Operator(std::string prot_selected) : DE_Operator() {
  init_files(prot_selection[prot_selected]);
  init_setup();
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
  scorefxn = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score1").c_str());
  std::cout << "start score: " << (*scorefxn)(*pose_) << std::endl;
  de_len = 10;
  start_res = 2;
  best_score = 1000000;
  init_two_stages_mover();

  // COMPLETE ABINITIO MOVER AT THE BOTTOM OF THE FILE
}

CalculateRmsdDistancePopulationPtr
DE_Operator::use_distances_strategy(std::string option) {
  switch (distances_map[option]) {
  case rmsd: {
    calculate_distances_popul = CalculateRmsdDistancePopulationPtr( new CalculateRmsdDistancePopulation(native_pose_, ffxn, ss, scorefxn, fit_radius));
    break;
  }
  case euclidean: {
    calculate_distances_popul = CalculateRmsdDistancePopulationPtr( new CalculateEuclideanDistancePopulation(native_pose_, ffxn, ss, scorefxn, fit_radius));
    break;
  }
  case euclidean_loop: {
    calculate_distances_popul = CalculateRmsdDistancePopulationPtr( new CalculateEuclideanLoopDistancePopulation(native_pose_, ffxn, ss, scorefxn, fit_radius));
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

  frag_mover = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::greedy_search, frag_opt);
  ffxn = boost::shared_ptr<FitFunction>( new PoseFragmentFunction(pose_, scorefxn, ss, frag_mover));
  //ffxn = boost::shared_ptr<FitFunction>( new PoseToyFunction(pose_, scorefxn, ss, frag_mover));
  local_search_frag_mover = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::stage3mover, frag_opt);
  local_search = boost::shared_ptr<LocalSearchIndividualMover>(
							       new LocalSearchIndividualMover(pose_, scorefxn, ss, local_search_frag_mover));

  ffxn->name_ = std::string("Score Function for stage " + stage_name);
  print_best = PrintBestIndividualPtr( new PrintBestIndividual(pose_, ffxn, ss, scorefxn));

  //calculate_distances_popul = use_distances_strategy("rmsd");
  calculate_distances_popul = use_distances_strategy( app_options.get<std::string>("Protocol.distance_strategy") );

  init_popul = boost::shared_ptr<InitPopulation>(new TwoStagesInitPopulation(ffxn, pose_, two_stages_mover, ss));
  //    init_popul = boost::shared_ptr<InitPopulation>(new PosePopulation(ffxn, pose_, ss));

  boost::shared_ptr<MoverDE> de = init_differential_evolution_protocol();

  de->Gmax = gmax_per_stage[stage_name];
  de->use_print_class = true;
  de->print_best_ind = print_best;
  de->calculate_distances_popul = calculate_distances_popul;

  set_current_population(de);

  de->popul = current_population;

  return de;
}

void
DE_Operator::init_two_stages_mover() {
  core::scoring::ScoreFunctionOP score1 = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score1").c_str());
  two_stages_mover = boost::shared_ptr<InitStagesMover>( new InitStagesMover(score1, frag_set_, frag_set_large));
}

void
DE_Operator::set_current_population( boost::shared_ptr<MoverDE> de ) {
  if (current_population.size() == 0) {
    current_population = de->popul;
  } else {
    CalculateRmsdDistancePopulation::DistancesResult result =
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
  ConfigurationDE conf = ConfigurationDE(app_options.get<int>("DE.NP"), app_options.get<double>("DE.CR"), app_options.get<double>("DE.F"), app_options.get<double>("Extra.fitrad") );
  std::string protocol_name = app_options.get<std::string>("Protocol.name");
  boost::shared_ptr<MoverDE> de ;
  switch (protocol_name_map[protocol_name]) {
  case Shared: {
    de = boost::shared_ptr<MoverDE>(new SharedMoverDE(conf, ffxn, init_popul));
    break;
  }
  case HybridShared: {
    de = boost::shared_ptr<MoverDE>(new SharedHybridMoverDE(conf, ffxn, init_popul, local_search));
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
DE_Operator::run() {
  std::string stage_name;
  std::cout << "run " << std::endl;
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
