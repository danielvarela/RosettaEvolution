
#ifndef ALGORITHMBUILDER_H
#define ALGORITHMBUILDER_H

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
#include "../Extra/Utils.hpp"
#include "OptionsMapInitializer.hpp"
#include "PopulStrategyInitializer.hpp"
#include "FitnessFunctionInitializer.hpp"

class AlgorithmBuilder
{
public:
  boost::property_tree::ptree app_options;
  boost::shared_ptr<FitFunction> ffxn;
  std::vector<Individual> current_population;
  boost::shared_ptr<LocalSearchIndividualMover> local_search;
  FragInsertionStrategy::FragOptions frag_opt;
  std::map<std::string, OptionsMapInitializer::protocol_name_enum> protocol_name_map;
  boost::shared_ptr<ConfigurationDE> conf;
  std::string protocol_name;


  AlgorithmBuilder(boost::property_tree::ptree app_options_in, boost::shared_ptr<FitFunction> ffxn_in,  std::vector<Individual> current_population_in,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in,  FragInsertionStrategy::FragOptions frag_opt_in);

  boost::shared_ptr<MoverDE> build(std::string default_protocol = "");

};


#endif /* ALGORITHMBUILDER_H */
