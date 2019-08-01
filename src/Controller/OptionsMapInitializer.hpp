#ifndef OPTIONSMAPINITIALIZER_H
#define OPTIONSMAPINITIALIZER_H

#include <string>
#include <map>
#include <vector>
#include <boost/tokenizer.hpp>

#include "../Extra/Utils.hpp"
#include "../Algorithm/DifferentialEvolutionMover.hpp"
#include "../Movers/PoseFunction.hpp"

namespace OptionsMapInitializer
{

  enum distances_enum {
    dist_error, rmsd, partial_rmsd, rmsd_native_diff, euclidean, euclidean_loop, euclidean_diff_abs, euclidean_partial_mario, euclidean_mario, mario_first_last
  };

  static std::map<std::string, distances_enum> distances_map_build() {
    std::map<std::string, distances_enum> distances_map;
    distances_map["error"] = dist_error;
    distances_map["rmsd"] = rmsd;
    distances_map["partial_rmsd"] = partial_rmsd;
    distances_map["rmsd_native_diff"] = rmsd_native_diff;
    distances_map["euclidean"] = euclidean;
    distances_map["euclidean_loop"] = euclidean_loop;
    distances_map["euclidean_diff_abs"] = euclidean_diff_abs;
    distances_map["euclidean_mario"] = euclidean_mario;
    distances_map["euclidean_partial_mario"] = euclidean_partial_mario;
    distances_map["mario_first_last"] = mario_first_last;
    return distances_map;
  }


  enum protocol_name_enum{
    protocol_error, Shared, HybridMover, HybridShared, ResetOldIndsHybrid ,  ResetOldIndsCrowdingHybrid , MPISeedsDE, MPIfasterCrowdingDE, MPIResetOldCrowdingDE, MPICrowdingDE, CrowdingDE, MPIMoverDE, HybridCrowdingDE
  };

  static std::map<std::string, protocol_name_enum> protocol_name_map_build(){
    std::map<std::string, protocol_name_enum> protocol_name_map;
    protocol_name_map["error"] = protocol_error;
    protocol_name_map["Shared"] = Shared;
    protocol_name_map["HybridShared"] = HybridShared;
    protocol_name_map["HybridMover"] = HybridMover;
    protocol_name_map["MPIMoverDE"] = MPIMoverDE;
    protocol_name_map["CrowdingDE"] = CrowdingDE;
    protocol_name_map["MPICrowdingDE"] = MPICrowdingDE;
    protocol_name_map["MPIfasterCrowdingDE"] =  MPIfasterCrowdingDE;
    protocol_name_map["MPISeedsDE"] = MPISeedsDE;
    protocol_name_map["MPIResetOldCrowdingDE"] = MPIResetOldCrowdingDE;
    protocol_name_map["HybridCrowdingDE"] = HybridCrowdingDE;
    protocol_name_map["ResetOld"] = ResetOldIndsHybrid;
    protocol_name_map["ResetOldCDE"] = ResetOldIndsCrowdingHybrid;
    return protocol_name_map;
  }

  enum init_popul_strategy_enum{
    popul_error, total_random, total_random_pose_based, random_pose_based, init_popul_with_stage
  };

  static std::map<std::string, init_popul_strategy_enum> init_popul_strategy_map_build() {
    std::map<std::string, init_popul_strategy_enum> init_popul_strategy_map;
    init_popul_strategy_map["error"] = popul_error;
    init_popul_strategy_map["total_random"] = total_random;
    init_popul_strategy_map["random_pose_based"] = random_pose_based;
    init_popul_strategy_map["total_random_pose_based"] = total_random_pose_based;
    init_popul_strategy_map["init_popul_with_stage"] = init_popul_with_stage;
    return init_popul_strategy_map;
  }

  static std::map<std::string, std::string> score_per_stage_build() {
    std::map<std::string, std::string> score_per_stage;

    score_per_stage["stage1"] = "score0";
    score_per_stage["stage2"] = "score1";
    score_per_stage["stage3"] = "score2";
    score_per_stage["stage4"] = "score3";
    //score_per_stage["stage4"] = "score4_smooth_cart";


    return score_per_stage;
  }

  static std::map<std::string, int> gmax_per_stage_build(const std::string& option = "default") {
    std::map<std::string, int> gmax_per_stage;
    if (option == "default") {
      gmax_per_stage["stage1"] = 30;
      gmax_per_stage["stage2"] = 100;
      gmax_per_stage["stage3"] = 100;
      gmax_per_stage["stage4"] = 100;
    } else {
      if (option == "short_test") {
	      gmax_per_stage["stage1"] = 10;
	      gmax_per_stage["stage2"] = 10;
	      gmax_per_stage["stage3"] = 10;
	      gmax_per_stage["stage4"] = 10;
      }
    }
    return gmax_per_stage;
  }



  static std::vector<std::string> parse_stages(const std::string& option) {
    std::string input_stages = option;
    trim(input_stages);
    std::vector<std::string> vec_stages;
    boost::char_separator<char> sep(",");
    boost::tokenizer<boost::char_separator<char> > tokens(input_stages, sep);
    for (boost::tokenizer<boost::char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it) {
    vec_stages.push_back(*it);
    }
    return vec_stages;
  }


  class CreateScoreFunctionStrategy
  {
  public:
    virtual core::scoring::ScoreFunctionOP build( const std::string& option) = 0;
    virtual CreateScoreFunctionStrategy* clone() const = 0;
    static boost::shared_ptr<OptionsMapInitializer::CreateScoreFunctionStrategy> get();
  };

  class RosettaDefaultFunctionCreator : public OptionsMapInitializer::CreateScoreFunctionStrategy
  {
  public:
    core::scoring::ScoreFunctionOP build(const std::string& option ) override;
    virtual RosettaDefaultFunctionCreator* clone() const {
      return new RosettaDefaultFunctionCreator(*this);
    }

  };

  class RosettaDensityFunctionCreator : public OptionsMapInitializer::CreateScoreFunctionStrategy
  {
  public:
     core::scoring::ScoreFunctionOP build(const std::string& option) override;
    virtual RosettaDensityFunctionCreator* clone() const {
      return new RosettaDensityFunctionCreator(*this);
    }
  };



};

#endif /* OPTIONSMAPINITIALIZER_H */
