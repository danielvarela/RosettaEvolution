
#ifndef POPULSTRATEGYINITIALIZER_H
#define POPULSTRATEGYINITIALIZER_H

#include <string>
#include <map>
#include <vector>
#include <boost/tokenizer.hpp>

#include "../Extra/Utils.hpp"
#include "../Algorithm/DifferentialEvolutionMover.hpp"
#include "../Movers/PoseFunction.hpp"
#include "../Movers/InitPopulation.hpp"
#include "../Movers/SecondStage.hpp"
#include "../Movers/FragInsertionMover.hpp"
#include "../Movers/CalculateRmsdDistancePopulation.hpp"
#include "../Movers/LocalSearchIndividualMover.hpp"
#include "FitnessFunctionInitializer.hpp"
#include "OptionsMapInitializer.hpp"

namespace OptionsMapInitializer
{

  class PopulStrategyInitializer
  {
  public:
    boost::shared_ptr<InitStagesMover> two_stages_mover;
    std::string init_popul_option;
    FragInsertionStrategy::FragOptions frag_opt;
    boost::shared_ptr<FitFunction> ffxn;
    double pose_disturb_degrees;
    std::map<std::string, OptionsMapInitializer::init_popul_strategy_enum> init_popul_strategy_map;

    PopulStrategyInitializer( boost::property_tree::ptree app_options, FragInsertionStrategy::FragOptions frag_opt_in);

    boost::shared_ptr<InitPopulation> build();
  };

}
#endif /* POPULSTRATEGYINITIALIZER_H */
