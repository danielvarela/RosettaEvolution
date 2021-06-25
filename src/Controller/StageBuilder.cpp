
#include "StageBuilder.hpp"

StageBuilder::StageBuilder(boost::property_tree::ptree app_options_default, core::scoring::ScoreFunctionOP scorefxn_default, FragInsertionStrategy::FragOptions frag_opt_default) {
  app_options = app_options_default;
  scorefxn = scorefxn_default;
  frag_opt = frag_opt_default;
  fit_radius = app_options.get<double>("Extra.fitrad");
  calculate_distances_option = app_options.get<std::string>("Protocol.distance_strategy") ;
  score_per_stage = OptionsMapInitializer::score_per_stage_build();
  gmax_per_stage = OptionsMapInitializer::gmax_per_stage_build();
}

bool
StageBuilder::initializer_frag_mover() {
  bool fragment_at_trials = false;
  boost::property_tree::ptree::const_assoc_iterator it = app_options.find("Fragments.strategy_at_trials");
  if(app_options.not_found() == it) {
    fragment_at_trials = true;
    std::string input_option = app_options.get<std::string>("Fragments.strategy_at_trials");
    frag_mover = frag_inserter_builder.get(input_option, frag_opt);
  } else {
    fragment_at_trials = false;
    std::string input_option = std::string("my_frag_insertion");
    frag_mover = frag_inserter_builder.get(input_option, frag_opt);
  }
  return fragment_at_trials;
}

bool
StageBuilder::initializer_frag_popul() {
  bool fragment_at_popul = false;
  boost::property_tree::ptree::const_assoc_iterator it = app_options.find("Fragments.strategy_at_population");
  if(app_options.not_found() == it) {
    fragment_at_popul = true;
    std::string input_option = app_options.get<std::string>("Fragments.strategy_at_population");
    popul_frag_mover = frag_inserter_builder.get(input_option, frag_opt);
  } else {
    //default option
    popul_frag_mover = frag_inserter_builder.get(std::string("stage_rosetta_mover"), frag_opt);
  }
  local_search = boost::shared_ptr<LocalSearchIndividualMover>(new LocalSearchIndividualMover(frag_opt.native_model, frag_opt.scorefxn, frag_opt.ss, popul_frag_mover));
  return fragment_at_popul;
}

PrintBestIndividualPtr
StageBuilder::initializer_print_best() {
  return PrintBestIndividualPtr( new PrintBestIndividual(frag_opt.native_model, ffxn, frag_opt.ss, frag_opt.scorefxn));
}

void
StageBuilder::update_current_population_to_stage_score( std::vector<Individual>& current_population ) {
  if (current_population.size() == 0) {
    //    current_population = de->popul;
    std::cout << "error, the current population was not initialize" << std::endl;
    exit(1);
  } else {
    //CalculateDistancePopulation::DistancesResult result =
    //  calculate_distances_popul->run(current_population);
    //std::vector<double> current_rmsd = result.rmsd_to_native;
    // rescore population

    boost::shared_ptr<MasterRosettaCalculator> mpi_calculator;
#if(MPI_ENABLED)
  
    ffxn = OptionsMapInitializer::FitnessFunctionInitializer(frag_opt, frag_mover, true).set_up();
    //mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterRosettaCalculator(ffxn));
    mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterScatterGather(ffxn, calculate_distances_popul));

  mpi_calculator->stage_name = ffxn->name();
#else
    ffxn = OptionsMapInitializer::FitnessFunctionInitializer(frag_opt, frag_mover, false).set_up();
#endif

#if(MPI_ENABLED)
    int frags_at_popul = 1;
    current_population = mpi_calculator->run(current_population, frags_at_popul);
#else
    for (int i = 0; i < current_population.size(); i++) {
      ffxn->score(current_population[i]);
    }
#endif
    std::cout << "end score popul to current score " << std::endl;
  }
}

void
StageBuilder::initializer_distances_calculator() {
  calculate_distances_popul = CalculateRmsdDistancePopulation(frag_opt.native_model, ffxn, frag_opt.ss, frag_opt.scorefxn, fit_radius).use_distances_strategy(calculate_distances_option);
  calculate_partial_rmsd = CalculateRmsdDistancePopulation(frag_opt.native_model, ffxn, frag_opt.ss, frag_opt.scorefxn, fit_radius).use_distances_strategy("partial_rmsd");
}

boost::shared_ptr<MoverDE>
StageBuilder::prepare_DE_for_stage(std::string stage_name, std::vector<Individual>& current_population) {
  scorefxn = core::scoring::ScoreFunctionFactory::create_score_function(std::string(score_per_stage[stage_name]).c_str());
  frag_opt.scorefxn = scorefxn;
  frag_opt.stage_name = stage_name;
  bool fragment_at_trials = initializer_frag_mover();
  ffxn = OptionsMapInitializer::FitnessFunctionInitializer(frag_opt, frag_mover, fragment_at_trials ).set_up();
  bool fragment_at_popul = initializer_frag_popul();
  ffxn->name_ = std::string("Score Function for stage " + stage_name);
  print_best = initializer_print_best();
  initializer_distances_calculator();
  //update_current_population_to_stage_score(current_population);
  boost::shared_ptr<MoverDE> de = AlgorithmBuilder(app_options, ffxn,  current_population,   local_search, frag_opt).build();
  de->Gmax = gmax_per_stage[stage_name];
  de->use_print_class = true;
  de->print_best_ind = print_best;
  de->calculate_distances_popul = calculate_distances_popul;
  de->calculator_partial_rmsd = boost::dynamic_pointer_cast<PartialRMSDcalculator>(calculate_partial_rmsd);
  de->popul = current_population;
  return de;
}
