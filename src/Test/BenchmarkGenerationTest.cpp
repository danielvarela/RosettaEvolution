
#include "BenchmarkGenerationTest.hpp"
#include <vector>

BenchmarkGenerationTest::BenchmarkGenerationTest(boost::shared_ptr<DE_Operator> app_operator_in) {
  app_operator = app_operator_in;
  app_operator->init_defaults();
  init_operator_for_test();
}

std::vector<Individual>
BenchmarkGenerationTest::start_popul() {
  std::vector<Individual> local_population(0);

  int NP = app_operator->app_options.get<int>("DE.NP");
  app_operator->init_popul = app_operator->initialize_init_popul_strategy("random_pose_based");

  app_operator->init_popul->apply(local_population, NP, app_operator->ffxn->D() );
  //  boost::shared_ptr<MoverDE> de = app_operator->init_differential_evolution_protocol("MPICrowdingDE");
  // std::cout << "Init Population " << std::endl;
  // for (int i = 0; i < NP; i++) {
  //   std::cout << "ind " << i << " : " << local_population[i].score << std::endl;
  // }
  return local_population;
}

void
BenchmarkGenerationTest::init_operator_for_test() {
  current_population = start_popul();
  app_operator->gmax_per_stage = OptionsMapInitializer::gmax_per_stage_build("short_test");
  StageBuilder stage_builder(app_operator->app_options,app_operator->scorefxn, app_operator->frag_opt);
  //  boost::shared_ptr<ConfigurationDE> conf = boost::shared_ptr<ConfigurationDE>( new ConfigurationDE(NP, 1.0, 0.025, NP, app_options.get<std::string>("Protocol.prot") ) );
  de = stage_builder.prepare_DE_for_stage( "stage4", current_population);
}

bool
BenchmarkGenerationTest::run() {
  timestamp_t t0 = get_timestamp();
  de->Gmax = 1;
  de->apply();
  current_population = de->popul;
  timestamp_t t1 = get_timestamp();
  double secs = (t1 - t0) / 1000000.0L;
  std::cout << "TOTAL BENCHMARK time : " << secs << std::endl;
  return true;
}
