
#include "OptionsMapInitializer.hpp"

boost::shared_ptr<OptionsMapInitializer::CreateScoreFunctionStrategy>
OptionsMapInitializer::CreateScoreFunctionStrategy::get() {
#if(USE_CRYO_EM)
  return  boost::shared_ptr<OptionsMapInitializer::CreateScoreFunctionStrategy>(  new OptionsMapInitializer::RosettaDensityFunctionCreator() );
#else
  return boost::shared_ptr<OptionsMapInitializer::CreateScoreFunctionStrategy>(   new OptionsMapInitializer::RosettaDefaultFunctionCreator() );
#endif
}

core::scoring::ScoreFunctionOP
OptionsMapInitializer::RosettaDefaultFunctionCreator::build(const std::string& option ) {
  return core::scoring::ScoreFunctionFactory::create_score_function(std::string(option).c_str());
}

core::scoring::ScoreFunctionOP
OptionsMapInitializer::RosettaDensityFunctionCreator::build(const std::string& option) {
  core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function(std::string(option).c_str());
  scorefxn->set_weight( core::scoring::elec_dens_fast, 30.0 );
  return scorefxn;
}
