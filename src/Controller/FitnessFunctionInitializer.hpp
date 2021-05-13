
#ifndef FITNESSFUNCTIONINITIALIZER_H
#define FITNESSFUNCTIONINITIALIZER_H

#include <boost/tokenizer.hpp>

#include "../Extra/Utils.hpp"
#include "../Algorithm/DifferentialEvolutionMover.hpp"
#include "../Movers/PoseFunction.hpp"
#include "../Movers/InitPopulation.hpp"
#include "../Movers/SecondStage.hpp"
#include "../Movers/FragInsertionMover.hpp"
#include "../Movers/CalculateRmsdDistancePopulation.hpp"
#include "../Movers/LocalSearchIndividualMover.hpp"
#include "../Extra/Utils.hpp"


namespace OptionsMapInitializer {
class FitnessFunctionInitializer
{
public:
FragInsertionStrategy::FragOptions frag_opt;
FragInsertionMoverPtr frag_mover;
bool fragment_at_trials;

FitnessFunctionInitializer(FragInsertionStrategy::FragOptions frag_opts_in, FragInsertionMoverPtr frag_mover_in, bool fragment_at_trials_in ) ;
boost::shared_ptr<FitFunction> set_up();
};

}


#endif /* FITNESSFUNCTIONINITIALIZER_H */
