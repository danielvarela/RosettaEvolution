
#ifndef LOCALSEARCHINDIVIDUALMOVER_H
#define LOCALSEARCHINDIVIDUALMOVER_H

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//Density
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/fragment/Frame.hh>

#include "DE_types.hpp"
#include "FitFunction.hpp"
#include "FragInsertionMover.hpp"
#include "PoseFunction.hpp"


class LocalSearchIndividualMover : public PoseFragmentFunction
{
public:

  LocalSearchIndividualMover(core::pose::PoseOP p, core::scoring::ScoreFunctionOP sfxn, std::string ss_in, FragInsertionMoverPtr frag_mover_) :   PoseFragmentFunction(p, sfxn, ss_in, frag_mover_) {
    frag_mover = frag_mover_;
  }

  double score(Individual& ind) {
    return apply(ind);
  }

  void reset_stats() {
    stats.clear();
  }

  std::string print_stats() override {


    int total_improves = 0;
    if (stats.find("improved_after_LS") != stats.end()) {
      total_improves = stats["improved_after_LS"];
    } else {
      total_improves = 0;
    }
    std::string output_string = "[LS_IMPROVED] " + std::to_string(total_improves);
    if (total_improves > 0) {
      output_string += " amt_improves " + std::to_string( stats["improved_amt_LS"] / stats["improved_after_LS"] );
      } else {
      output_string += " amt_improves 0";
    }

    if (stats.find("total_tries_LS") != stats.end()) {
      output_string += " total_tries " + std::to_string(stats["total_tries_LS"]);
    }
    if (stats.find("accepted_fragments_LS") != stats.end()) {
      output_string += " accepted_fragments " + std::to_string(stats["accepted_fragments_LS"]);
    }


    return output_string;
  }

  virtual double apply(Individual& ind) {
    double score_at_individual = ind.score;
    std::vector<double> init_vars = ind.vars;
    fill_pose(pose_, ind, ss);
    double result = SCORE_ERROR_FIXED + (*scorefxn)(*pose_);
    pose_to_ind(pose_, ind);
    ind.score = result;
    double result_after_frags = run_frag_mover(*pose_);
    if (result_after_frags < result) {
      //from pose to ind
      pose_to_ind(pose_, ind);
      ind.score = result_after_frags;
      result = result_after_frags;
      if (stats.find("improved_after_LS") != stats.end()) {
	stats["improved_after_LS"]++;

	stats["improved_amt_LS"] += ( result_after_frags - result );
      } else {
	stats["improved_after_LS"] = 1;
	stats["improved_amt_LS"] = 0;
      }
    }

    if (stats.find("total_tries") != stats.end()) {
      if (stats.find("total_tries_LS") == stats.end()) {
	stats["total_tries_LS"] = 0;
      }
      stats["total_tries_LS"] += stats["total_tries"];
      stats.erase("total_tries");
    }

    if (stats.find("accepted_fragments") != stats.end()) {
      stats["accepted_fragments_LS"] = stats["accepted_fragments"];
      stats.erase("accepted_fragments");
    }

    return result;
  }
};


#endif /* LOCALSEARCHINDIVIDUALMOVER_H */
