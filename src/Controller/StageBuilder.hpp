
#ifndef STAGEBUILDER_H
#define STAGEBUILDER_H

#include <string>
#include <map>
#include <vector>
#include <boost/tokenizer.hpp>
#include "../Extra/Utils.hpp"
#include "OptionsMapInitializer.hpp"
#include "AlgorithmBuilder.hpp"
#include "FragmentInserterBuilder.hpp"

class StageBuilder
{
public:
  boost::property_tree::ptree app_options;
  core::scoring::ScoreFunctionOP scorefxn;
  FragInsertionStrategy::FragOptions frag_opt;
  CalculateDistancePopulationPtr calculate_distances_popul;

private:
  FragInsertionMoverPtr frag_mover;
  FragInsertionMoverPtr popul_frag_mover;
  FragmentInserterBuilder frag_inserter_builder;
  boost::shared_ptr<FitFunction> ffxn;
  boost::shared_ptr<LocalSearchIndividualMover> local_search;
  PrintBestIndividualPtr print_best;
 CalculateDistancePopulationPtr calculate_partial_rmsd;
  double fit_radius;
  std::string calculate_distances_option;
  std::map<std::string, std::string> score_per_stage;
  std::map<std::string, int> gmax_per_stage;
//METHODS
public:
  StageBuilder(boost::property_tree::ptree app_options_default, core::scoring::ScoreFunctionOP scorefxn_default, FragInsertionStrategy::FragOptions frag_opt_default);

  bool initializer_frag_mover();
  bool initializer_frag_popul();
  PrintBestIndividualPtr initializer_print_best() ;
  void update_current_population_to_stage_score( std::vector<Individual>& current_population ) ;
  void initializer_distances_calculator();

  boost::shared_ptr<MoverDE> prepare_DE_for_stage(std::string stage_name, std::vector<Individual>& current_population);

};



#endif /* STAGEBUILDER_H */
