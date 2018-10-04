
#ifndef PRINTBESTINDIVIDUAL_H
#define PRINTBESTINDIVIDUAL_H

#include "DE_types.hpp"

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>

#include "FitFunction.hpp"
#include "PoseFunction.hpp"

class PrintBestIndividual {
public:
  core::pose::PoseOP pose_;
  int cnt;
  int uniq_id;
  FitFunctionPtr scorefxn;
  std::string ss;
  core::pose::PoseOP native_;
  core::scoring::ScoreFunctionOP pose_score;
  core::pose::PoseOP previous_best;

  PrintBestIndividual(const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore);

  std::string
  print(Individual ind);

  std::string
  print(Individual ind, std::string id);

  std::string
  print_energies(Individual ind);

};

typedef boost::shared_ptr<PrintBestIndividual> PrintBestIndividualPtr;



#endif /* PRINTBESTINDIVIDUAL_H */
