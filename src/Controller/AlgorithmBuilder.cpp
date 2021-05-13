
#include "AlgorithmBuilder.hpp"

AlgorithmBuilder::AlgorithmBuilder(boost::property_tree::ptree app_options_in, boost::shared_ptr<FitFunction> ffxn_in,  std::vector<Individual> current_population_in,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in,  FragInsertionStrategy::FragOptions frag_opt_in) {
  app_options = app_options_in;
  ffxn = ffxn_in;
  current_population = current_population_in;
  local_search = local_search_in;
  frag_opt = frag_opt_in;
  protocol_name_map = OptionsMapInitializer::protocol_name_map_build();
  conf = boost::shared_ptr<ConfigurationDE>( new ConfigurationDE(app_options.get<int>("DE.NP"), app_options.get<double>("DE.CR"), app_options.get<double>("DE.F"), app_options.get<double>("Extra.fitrad") , app_options.get<std::string>("Protocol.prot")) );
}

boost::shared_ptr<MoverDE>
AlgorithmBuilder::build(std::string default_protocol) {
  //     boost::shared_ptr<MoverDE> de = boost::shared_ptr<MoverDE>(new SharedMoverDE(ffxn, init_popul));
  protocol_name = default_protocol;
  if (protocol_name.size() <= 1) {
    protocol_name = app_options.get<std::string>("Protocol.name");
  }

  std::string crossover_strategy = app_options.get<std::string>("Protocol.crossover_strategy");
  std::string select_parents_strategy = app_options.get<std::string>("Protocol.select_parents_strategy");

  boost::shared_ptr<MoverDE> de ;
  switch (protocol_name_map[protocol_name]) {
  case OptionsMapInitializer::protocol_name_enum::Shared: {
    de = boost::shared_ptr<MoverDE>(new SharedMoverDE(app_options, ffxn, current_population));
    break;
  }
  case OptionsMapInitializer::protocol_name_enum::HybridMover: {
    de = boost::shared_ptr<MoverDE>(new HybridMoverDE(app_options, ffxn, current_population, local_search));
    break;
  }
  case OptionsMapInitializer::protocol_name_enum::HybridShared: {
    de = boost::shared_ptr<MoverDE>(new SharedHybridMoverDE(app_options, ffxn, current_population, local_search));
    break;
  }
  case OptionsMapInitializer::protocol_name_enum::CrowdingDE: {
    de = boost::shared_ptr<MoverDE>(new CrowdingMoverDE(app_options, ffxn, current_population));
    de->only_de = false;
    break;
  }
  case OptionsMapInitializer::protocol_name_enum::MPIMoverDE: {
    de = boost::shared_ptr<MoverDE>(new CrowdingMoverDE(app_options, ffxn, current_population));
    de->only_de = true;
    break;
  }
  case OptionsMapInitializer::protocol_name_enum::MPISeedsDE: {
    de = boost::shared_ptr<MoverDE>(new MPISeedsSantosMoverDE(app_options, ffxn, current_population, local_search));
    break;
  }
  case OptionsMapInitializer::protocol_name_enum::MPICrowdingDE: {
    de = boost::shared_ptr<MoverDE>(new MPICrowdingMoverDE(app_options, ffxn, current_population, local_search));
    break;
  }
  case OptionsMapInitializer::protocol_name_enum::MPIfasterCrowdingDE: {
    de = boost::shared_ptr<MoverDE>(new MPIfastCrowdingMoverDE(app_options, ffxn, current_population, local_search));
    break;
  }
  case OptionsMapInitializer::protocol_name_enum::HybridCrowdingDE: {
    de = boost::shared_ptr<MoverDE>(new CrowdingHybridMoverDE(app_options, ffxn, current_population, local_search));
    break;
  }
  default:
    std::cout << "error with Protocol Name" << std::endl;
    exit(1);
    break;
  }

  de->init_popul_ = OptionsMapInitializer::PopulStrategyInitializer(app_options, frag_opt).build();
  de->select_parents_strategy = select_parents_strategy;
  de->crossover_strategy = crossover_strategy;
  return de;

}
