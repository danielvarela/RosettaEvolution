
#include "CryoProtocolTest.hpp"
#include "../Controller/StageBuilder.hpp"

#include <core/scoring/Energies.hh>

#include <protocols/hybridization/BackboneTorsionPerturbation.hh>

CryoProtocolTest::CryoProtocolTest(boost::shared_ptr<DE_Operator> app_operator_in) {
  app_operator = app_operator_in;
  app_operator->init_defaults();
  init_operator_for_test();
}

std::vector<Individual>
CryoProtocolTest::start_popul() {
  std::vector<Individual> local_population(0);
  int NP = app_operator->app_options.get<int>("DE.NP");
  app_operator->init_popul = app_operator->initialize_init_popul_strategy("random_pose_based");
  app_operator->init_popul->apply(local_population, NP, app_operator->ffxn->D() );
  app_operator->gmax_per_stage = OptionsMapInitializer::gmax_per_stage_build("short_test");
  StageBuilder stage_builder(app_operator->app_options,app_operator->scorefxn, app_operator->frag_opt);
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
CryoProtocolTest::init_operator_for_test() {
  using namespace core::scoring;
  current_population = start_popul();
  app_operator->gmax_per_stage = OptionsMapInitializer::gmax_per_stage_build("short_test");
  StageBuilder stage_builder(app_operator->app_options,app_operator->scorefxn, app_operator->frag_opt);
  ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart");
  scorefxn->set_weight( core::scoring::elec_dens_fast, 30.0 );
  app_operator->scorefxn = scorefxn;
  std::cout << "start prepare popul" << std::endl;
  de = stage_builder.prepare_DE_for_stage( "stage4", current_population);
  current_population = de->popul;
  protocols::electron_density::SetupForDensityScoringMoverOP dockindens = protocols::electron_density::SetupForDensityScoringMoverOP( new protocols::electron_density::SetupForDensityScoringMover );
  core::pose::PoseOP pose_ = app_operator->native_pose_->clone();
  std::cout << "current rest population " << std::endl;
  app_operator->scorefxn->set_weight( core::scoring::elec_dens_fast, 30.0 );
  pose_->energies().show_total_headers(std::cout);
  for (int i = 0; i < current_population.size(); i++) {
    dockindens->apply(*pose_);
    de->scfxn->fill_pose(pose_, current_population[i], app_operator->ss);

    double score_no_density = app_operator->scorefxn->score(*pose_);
    std::cout << "FIRST ind " << i << " : " << (-1 * SCORE_ERROR_FIXED) +  current_population[i].score << " ";
    std::cout << "without dendisty " << score_no_density << " " << std::endl;
    pose_->energies().show_totals(std::cout);
    std::cout << std::endl;
    core::scoring::EnergyMap energy_map = pose_->energies().total_energies();
    std::cout << "density " << std::to_string(energy_map[core::scoring::ScoreType::elec_dens_fast]) << " " << std::endl;
    //    pose_->dump_pdb("pose_"+std::to_string(i) + "_.pdb");
  }
  scfxn = de->scfxn;
}

class RosettaCartesianHybridize
{
public:
  RosettaCartesianHybridize() {
	// CartesianHybridize(
	// 	utility::vector1 < core::pose::PoseOP > const & templates_in,
	// 	utility::vector1 < core::Real > const & template_wts_in,
	// 	utility::vector1 < protocols::loops::Loops > const & template_chunks_in,
	// 	utility::vector1 < protocols::loops::Loops > const & template_contigs_in,
	// 	core::fragment::FragSetOP fragments9_in );

    using namespace protocols::hybridization;

     
  }
  
  virtual ~RosettaCartesianHybridize();
};


bool
CryoProtocolTest::run() {
  using namespace protocols::hybridization;
  using namespace core::scoring;
  
  protocols::electron_density::SetupForDensityScoringMoverOP dockindens = protocols::electron_density::SetupForDensityScoringMoverOP( new protocols::electron_density::SetupForDensityScoringMover );
  core::pose::PoseOP pose_ = app_operator->native_pose_->clone();
  dockindens->apply(*pose_);

  BackboneTorsionPerturbation bb_torsion;
  bb_torsion.set_scorefunction(app_operator->scorefxn);

  CartesianHybridize cart_hybrid;

  ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score4_smooth_cart");
  scorefxn->set_weight( core::scoring::elec_dens_fast, 30.0 );
  
  cart_hybrid.set_scorefunction(scorefxn);
  cart_hybrid.set_no_global_frame(false);
  cart_hybrid.set_linmin_only(false);
  cart_hybrid.set_cenrot(false);
  cart_hybrid.set_cartfrag_overlap(7);
  cart_hybrid.set_temperature(2.0);
  cart_hybrid.set_max_insertion(50);

  
  std::cout << "start score " << app_operator->scorefxn->score(*pose_) << std::endl;
  cart_hybrid.apply(*pose_); 
  //bb_torsion.apply(*pose_);
  //bb_torsion.perturb(*pose_, 2.0);
  std::cout << "smooth_cart score " << scorefxn->score(*pose_) << std::endl;
  std::cout << "end score " << app_operator->scorefxn->score(*pose_) << std::endl;

  for (int i = 0; i < current_population.size(); i++) {
    dockindens->apply(*pose_);
    de->scfxn->fill_pose(pose_, current_population[i], app_operator->ss);


    std::cout << " POSE " << i << std::endl;
  std::cout << "start score " << app_operator->scorefxn->score(*pose_) << std::endl;
    cart_hybrid.apply(*pose_); 
    double score_no_density = app_operator->scorefxn->score(*pose_);

  }
  
  double start_score = current_population[0].score;

 
  //scfxn->score(current_population[0]);
  double final_score = current_population[0].score;
  if (start_score < final_score) {
    std::cout << "final score " << final_score << " > " << start_score << std::endl;
    return false;
  }
  return true;
}


