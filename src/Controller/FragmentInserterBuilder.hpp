#ifndef FRAGMENTINSERTERBUILDER_H
#define FRAGMENTINSERTERBUILDER_H


#include "../Extra/Utils.hpp"
#include "../Movers/FragInsertionMover.hpp"
#include <boost/property_tree/ptree.hpp>

class FragmentInserterBuilder
{
public:
  FragmentInserterBuilder() {
    fragment_insertion_strategy_map["error"] = frag_error;
    fragment_insertion_strategy_map["my_frag_insertion"] = my_frag_insertion;
    fragment_insertion_strategy_map["greedy_search"] = greedy_search;
    fragment_insertion_strategy_map["no_greedy_search"] = no_greedy_search;
    fragment_insertion_strategy_map["hybrid_mover"] = hybrid_mover;
    fragment_insertion_strategy_map["stage_rosetta_mover"] = stage_rosetta_mover;
    fragment_insertion_strategy_map["ILS_as_julia"] = ILS_as_julia;
  }

  enum fragment_insertion_strategy_enum {
    frag_error, my_frag_insertion, no_greedy_search, greedy_search, hybrid_mover, stage_rosetta_mover, ILS_as_julia
  };
  std::map<std::string, fragment_insertion_strategy_enum> fragment_insertion_strategy_map;

  FragInsertionStrategy::FragOptions
  init_frag_files(Protinfo prot_info, boost::property_tree::ptree app_options  );

  void
  init_density_map(core::pose::Pose& ipose, std::string prot_name, int resolution );

  boost::shared_ptr<FragInsertionMover> get(std::string option, FragInsertionStrategy::FragOptions frag_opt);

};


#endif /* FRAGMENTINSERTERBUILDER_H */
