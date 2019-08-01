
#include "InvididualEvaluationTest.hpp"
#include "../Controller/StageBuilder.hpp"

#include <core/scoring/Energies.hh>
IndividualEvaluationTest::IndividualEvaluationTest(boost::shared_ptr<DE_Operator> app_operator_in) {
  app_operator = app_operator_in;
  app_operator->init_defaults();
  init_operator_for_test();
}

std::vector<Individual>
IndividualEvaluationTest::start_popul() {
  std::vector<Individual> local_population(0);
  int NP = app_operator->app_options.get<int>("DE.NP");
  app_operator->init_popul = app_operator->initialize_init_popul_strategy("random_pose_based");

  app_operator->init_popul->apply(local_population, NP, app_operator->ffxn->D() );
  //  boost::shared_ptr<MoverDE> de = app_operator->init_differential_evolution_protocol("MPICrowdingDE");
  app_operator->gmax_per_stage = OptionsMapInitializer::gmax_per_stage_build("short_test");
  StageBuilder stage_builder(app_operator->app_options,app_operator->scorefxn, app_operator->frag_opt);
  //  boost::shared_ptr<ConfigurationDE> conf = boost::shared_ptr<ConfigurationDE>( new ConfigurationDE(NP, 1.0, 0.025, NP, app_options.get<std::string>("Protocol.prot") ) );
  protocols::electron_density::SetupForDensityScoringMoverOP dockindens = protocols::electron_density::SetupForDensityScoringMoverOP( new protocols::electron_density::SetupForDensityScoringMover );
  core::pose::PoseOP pose_ = app_operator->native_pose_->clone();
  std::cout << "current rest population " << std::endl;
  app_operator->scorefxn->set_weight( core::scoring::elec_dens_fast, 30.0 );
  std::cout << "Init Population " << std::endl;
  for (int i = 0; i < local_population.size(); i++) {
    app_operator->ffxn->fill_pose(pose_, local_population[i], app_operator->ss);
    dockindens->apply(*pose_);
    double score_no_density = app_operator->scorefxn->score(*pose_);
    local_population[i].score = (SCORE_ERROR_FIXED) +  score_no_density;
    std::cout << "ind " << i << " : " << local_population[i].score << std::endl;
  }
  return local_population;
}

void
IndividualEvaluationTest::init_operator_for_test() {
  current_population = start_popul();
  app_operator->gmax_per_stage = OptionsMapInitializer::gmax_per_stage_build("short_test");
  StageBuilder stage_builder(app_operator->app_options,app_operator->scorefxn, app_operator->frag_opt);
  //  boost::shared_ptr<ConfigurationDE> conf = boost::shared_ptr<ConfigurationDE>( new ConfigurationDE(NP, 1.0, 0.025, NP, app_options.get<std::string>("Protocol.prot") ) );
  std::cout << "start prepare popul" << std::endl;
  de = stage_builder.prepare_DE_for_stage( "stage4", current_population);
  current_population = de->popul;
  protocols::electron_density::SetupForDensityScoringMoverOP dockindens = protocols::electron_density::SetupForDensityScoringMoverOP( new protocols::electron_density::SetupForDensityScoringMover );
  core::pose::PoseOP pose_ = app_operator->native_pose_->clone();
  std::cout << "current rest population " << std::endl;
  app_operator->scorefxn->set_weight( core::scoring::elec_dens_fast, 30.0 );
  pose_->energies().show_total_headers(std::cout);
  for (int i = 0; i < current_population.size(); i++) {
    de->scfxn->fill_pose(pose_, current_population[i], app_operator->ss);
    dockindens->apply(*pose_);
    double score_no_density = app_operator->scorefxn->score(*pose_);
    std::cout << "FIRST ind " << i << " : " << (-1 * SCORE_ERROR_FIXED) +  current_population[i].score << " ";
    std::cout << "without dendisty " << score_no_density << " " << std::endl;
    pose_->energies().show_totals(std::cout);
    std::cout << std::endl;
    core::scoring::EnergyMap energy_map = pose_->energies().total_energies();
    std::cout << "density " << std::to_string(energy_map[core::scoring::ScoreType::elec_dens_fast]) << " " << std::endl;
    pose_->dump_pdb("pose_"+std::to_string(i) + "_.pdb");
  }

  scfxn = de->scfxn;
}

bool
IndividualEvaluationTest::run() {
  double start_score = current_population[0].score;
  scfxn->score(current_population[0]);
  double final_score = current_population[0].score;
  if (start_score < final_score) {
    std::cout << "final score " << final_score << " > " << start_score << std::endl;
    return false;
  }
  return true;
}
