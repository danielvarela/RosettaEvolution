#include "NeighborhoodSelectionTest.hpp"
#include <vector>

NeighborhoodSelectionTest::NeighborhoodSelectionTest(boost::shared_ptr<DE_Operator> app_operator_in) {
  app_operator = app_operator_in;
  app_operator->init_defaults();
  init_operator_for_test();
}

std::vector<Individual>
NeighborhoodSelectionTest::start_popul() {
  std::vector<Individual> local_population(0);
  int NP = app_operator->app_options.get<int>("DE.NP");
  app_operator->init_popul = app_operator->initialize_init_popul_strategy("random_pose_based");

  app_operator->init_popul->apply(local_population, NP, app_operator->ffxn->D() );
  //  boost::shared_ptr<MoverDE> de = app_operator->init_differential_evolution_protocol("MPICrowdingDE");
  std::cout << "Init Population " << std::endl;
  for (int i = 0; i < NP; i++) {
    std::cout << "ind " << i << " : " << local_population[i].score << std::endl;
  }
  return local_population;
}

void
NeighborhoodSelectionTest::init_operator_for_test() {
  current_population = start_popul();
  app_operator->gmax_per_stage = OptionsMapInitializer::gmax_per_stage_build("short_test");
  StageBuilder stage_builder(app_operator->app_options,app_operator->scorefxn, app_operator->frag_opt);
  //  boost::shared_ptr<ConfigurationDE> conf = boost::shared_ptr<ConfigurationDE>( new ConfigurationDE(NP, 1.0, 0.025, NP, app_options.get<std::string>("Protocol.prot") ) );
  de = stage_builder.prepare_DE_for_stage( "stage4", current_population);
  FitFunctionPtr ffxn = de->scfxn;
  boost::shared_ptr<PoseFragmentFunction> aux_ffxn = boost::dynamic_pointer_cast<PoseFragmentFunction>(ffxn);
  simple_ffxn = boost::shared_ptr<PoseScoreFunction>( new PoseScoreFunction(app_operator->pose_, app_operator->scorefxn, aux_ffxn->ss, aux_ffxn->frag_mover));
  for (int i = 0; i < current_population.size(); i++) {
    simple_ffxn->score(current_population[i]);
  }
  mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterRosettaCalculator(ffxn));


}

bool
NeighborhoodSelectionTest::run() {
  for (int i = 0; i < current_population.size(); i++) {
    int target_ind = de->calculate_distances_popul->find_nearest(current_population[i], current_population);
    if (target_ind != i) {
      return false;
    }
  }

  return true;
}
