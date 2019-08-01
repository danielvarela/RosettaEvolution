
#ifndef SECONDSTAGE_HPP
#define SECONDSTAGE_HPP

#include "../Algorithm/DE_types.hpp"
#include "FitFunction.hpp"
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pack/pack_rotamers.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <core/conformation/util.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/conformation/util.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/SingleResidueFragData.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>

#include <core/fragment/BBTorsionAndAnglesSRFD.hh>

class SecondStage
{
public:
  SecondStage();

  void apply(core::pose::PoseOP p);

protected:
  protocols::simple_moves::SwitchResidueTypeSetMoverOP to_fa;
  protocols::simple_moves::SwitchResidueTypeSetMoverOP to_cen;
  core::optimization::AtomTreeMinimizerOP rbminimizer;
  core::optimization::MinimizerOptionsOP  options_rb ;
  core::kinematics::MoveMapOP mm_rb;
  //  protocols::simple_moves::PackRotamersMoverOP	pack_mover ;
  core::scoring::ScoreFunctionOP densonly;
};

#endif // SECONDSTAGE_HPP
