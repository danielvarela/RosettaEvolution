
#include "FitnessFunctionInitializer.hpp"

OptionsMapInitializer::FitnessFunctionInitializer::FitnessFunctionInitializer(FragInsertionStrategy::FragOptions frag_opts_in, FragInsertionMoverPtr frag_mover_in, bool fragment_at_trials_in ) {
  frag_opt = frag_opts_in;
  frag_mover = frag_mover_in;
  fragment_at_trials = fragment_at_trials_in;
}

boost::shared_ptr<FitFunction>
OptionsMapInitializer::FitnessFunctionInitializer::set_up() {
  boost::shared_ptr<FitFunction> ffxn;
  if (fragment_at_trials) {
    ffxn = boost::shared_ptr<FitFunction>( new PoseFragmentFunction(frag_opt.native_model, frag_opt.scorefxn, frag_opt.ss, frag_mover));
#if(USE_CRYO_EM)
    std::cout << "ENTRA EN CRYO_EM" << std::endl;
    ffxn = boost::shared_ptr<FitFunction>( new PoseDensityFragmentFunction(frag_opt.native_model, frag_opt.scorefxn, frag_opt.ss, frag_mover));
#endif

    ffxn->name_ = frag_opt.stage_name;
  } else {
    //    ffxn = boost::shared_ptr<FitFunction>( new PoseToyFunction(pose_, scorefxn, ss, frag_mover));
    ffxn = boost::shared_ptr<FitFunction>( new PoseScoreFunction(frag_opt.native_model, frag_opt.scorefxn, frag_opt.ss, frag_mover));
#if(USE_CRYO_EM)
    ffxn = boost::shared_ptr<FitFunction>( new PoseDensityFunction(frag_opt.native_model, frag_opt.scorefxn, frag_opt.ss, frag_mover));
#endif
    ffxn->name_ = frag_opt.stage_name;
  }
  return ffxn;
}
