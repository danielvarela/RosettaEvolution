
#include "PopulStrategyInitializer.hpp"

OptionsMapInitializer::PopulStrategyInitializer::PopulStrategyInitializer( boost::property_tree::ptree app_options, FragInsertionStrategy::FragOptions frag_opt_in) {
  frag_opt = frag_opt_in;
  core::scoring::ScoreFunctionOP score1 = OptionsMapInitializer::CreateScoreFunctionStrategy::get()->build(std::string("score1"));
  two_stages_mover = boost::shared_ptr<InitStagesMover>( new InitStagesMover(score1, frag_opt.frag_set_, frag_opt.frag_set_large));
  frag_opt.scorefxn = score1;
  frag_opt.stage_name = "score1";
  ffxn = OptionsMapInitializer::FitnessFunctionInitializer(frag_opt, two_stages_mover, false).set_up();
  init_popul_option = app_options.get<std::string>("Protocol.init_strategy");
  pose_disturb_degrees = app_options.get<double>("Protocol.disturb_degrees");
  init_popul_strategy_map = OptionsMapInitializer::init_popul_strategy_map_build();
}

boost::shared_ptr<InitPopulation>
OptionsMapInitializer::PopulStrategyInitializer::build() {
  boost::shared_ptr<InitPopulation> init_popul_in_options_file;
  switch (init_popul_strategy_map[init_popul_option]) {
  case OptionsMapInitializer::init_popul_strategy_enum::total_random: {
    init_popul_in_options_file = boost::shared_ptr<InitPopulation>(new InitPopulation(ffxn));
    break;
  }
  case OptionsMapInitializer::init_popul_strategy_enum::total_random_pose_based: {
    init_popul_in_options_file = boost::shared_ptr<InitPopulation>(new PoseTotalRandomInitPopulation(ffxn, frag_opt.native_model, frag_opt.ss));
    break;
  }
  case OptionsMapInitializer::init_popul_strategy_enum::random_pose_based: {
    init_popul_in_options_file = boost::shared_ptr<InitPopulation>(new PosePopulation(ffxn, frag_opt.native_model, frag_opt.ss,   pose_disturb_degrees ));
    break;
  }
  case OptionsMapInitializer::init_popul_strategy_enum::init_popul_with_stage: {
    init_popul_in_options_file = boost::shared_ptr<InitPopulation>(new TwoStagesInitPopulation(ffxn, frag_opt.native_model, two_stages_mover, frag_opt.ss));
    break;
  }
  default:
    std::cout << "problem selecting the init population strategy" << std::endl;
    exit(1);
    break;
  }
  return init_popul_in_options_file;
}
