
#include "SuperDebugTest.hpp"
#include <vector>

SuperDebugTest::SuperDebugTest(boost::shared_ptr<DE_Operator> app_operator_in) {
  app_operator = app_operator_in;
  app_operator->init_defaults();
  init_operator_for_test();
}

std::vector<Individual>
SuperDebugTest::start_popul() {
  std::vector<Individual> local_population(0);
  int NP = app_operator->app_options.get<int>("DE.NP");
  app_operator->init_popul = app_operator->initialize_init_popul_strategy("random_pose_based");
  app_operator->init_popul->apply(local_population, NP, app_operator->ffxn->D() );
  //  boost::shared_ptr<MoverDE> de = app_operator->init_differential_evolution_protocol("MPICrowdingDE");
  return local_population;
}

void
SuperDebugTest::init_operator_for_test() {
  current_population = start_popul();
  app_operator->gmax_per_stage = OptionsMapInitializer::gmax_per_stage_build("short_test");
  StageBuilder stage_builder(app_operator->app_options,app_operator->scorefxn, app_operator->frag_opt);
  //  boost::shared_ptr<ConfigurationDE> conf = boost::shared_ptr<ConfigurationDE>( new ConfigurationDE(NP, 1.0, 0.025, NP, app_options.get<std::string>("Protocol.prot") ) );
  de = stage_builder.prepare_DE_for_stage( "stage4", current_population);
  FitFunctionPtr ffxn = de->scfxn;
  boost::shared_ptr<PoseFragmentFunction> aux_ffxn = boost::dynamic_pointer_cast<PoseFragmentFunction>(ffxn);
  simple_ffxn = boost::shared_ptr<PoseScoreFunction>( new PoseScoreFunction(app_operator->pose_, app_operator->scorefxn, aux_ffxn->ss, aux_ffxn->frag_mover));
  native_pose = app_operator->native_pose_;

  make_population_straight();
  insert_native_pose_at_popul();

  for (int i = 0; i < current_population.size(); i++) {
    simple_ffxn->score(current_population[i]);
  }

  std::cout << "Init Population " << std::endl;
  for (int i = 0; i < current_population.size(); i++) {
    std::cout << "ind " << i << " : " << current_population[i].score << std::endl;
  }


  mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterRosettaCalculator(ffxn));
}

void SuperDebugTest::insert_native_pose_at_popul() {
  simple_ffxn->pose_to_ind(native_pose, current_population[0]);
}


void
SuperDebugTest::make_population_straight() {
  for (int i = 0; i < current_population.size(); i++) {
    for (int j = 0; j < current_population[i].vars.size(); j++) {
      current_population[i].vars[j] = 1.0;
    }
  }

}

bool
SuperDebugTest::run() {
  std::vector<Individual> test_population = current_population;
  de->popul = current_population;
  de->Gmax = 1;
  de->apply();
  current_population = de->popul;

  std::cout << "[RANGE] " << std::endl;
  for (int i = 0; i < test_population.size(); i++) {
    if (i < test_population.size() - 1) {
      std::cout << (-1 * SCORE_ERROR_FIXED) +  current_population[i].score << " ";
    } else {
      std::cout << (-1 * SCORE_ERROR_FIXED) +  current_population[i].score << std::endl;
    }
  }


  for (int i = 1; i < test_population.size(); i++) {
    if (test_population[0].score < current_population[i].score) {
      std::cout << test_population[0].score << " < " <<  current_population[i].score << std::endl;
      //return false;
    }
  }

  return true;
}
