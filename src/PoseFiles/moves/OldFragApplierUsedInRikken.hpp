
#ifndef FRAGAPPLIER_HPP
#define FRAGAPPLIER_HPP

#include <protocols/hybridization/InsertChunkMover.hh>
#include <protocols/hybridization/CartesianSampler.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/kinematics_util.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>

#include <core/scoring/ScoreFunction.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <protocols/moves/Mover.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>
#include <numeric/random/WeightedSampler.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <boost/unordered/unordered_map.hpp>


using namespace core;
using namespace core::scoring;
using namespace core::scoring::constraints;

class FragApplier
{
public:
  FragApplier() {

    overlap_ = 2;
  }

  // transform fragment
  void
  apply_transform(
		  core::pose::Pose &frag,
		  core::Vector const &preT,
		  core::Vector const &postT,
		  numeric::xyzMatrix< core::Real > const &R
		  );


  // apply endpoint csts from pose to frag
  void
  apply_fragcsts(
		 core::pose::Pose &working_frag,
		 core::pose::Pose const &pose,
		 core::Size startpos
		 );


  // get frag->pose transform, return RMS
  // default overlap_ = 2
  core::Real
  get_transform(
		core::pose::Pose const &pose,
		core::pose::Pose const &frag,
		core::Size startpos,
		core::Vector &preT,
		core::Vector &postT,
		numeric::xyzMatrix< core::Real > &R
		);

  bool apply_frame (
		    core::pose::Pose & pose /*, core::fragment::Frame frame*/
		    );

protected:
  int overlap_;
};


#endif // FRAGAPPLIER_HPP
