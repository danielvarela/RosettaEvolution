
#ifndef FRAGINSERTIONMOVER_H
#define FRAGINSERTIONMOVER_H

#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MonteCarlo.hh>
#include "FitFunction.hpp"
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <protocols/hybridization/HybridizeProtocolCreator.hh>
#include <protocols/hybridization/HybridizeProtocol.hh>
#include <protocols/hybridization/FoldTreeHybridize.hh>
#include <protocols/hybridization/CartesianHybridize.hh>

#include <protocols/hybridization/CartesianSampler.hh>
#include <protocols/hybridization/CartesianSampler.fwd.hh>
#include <protocols/hybridization/TemplateHistory.hh>
#include <protocols/hybridization/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/hybridization/HybridizeProtocol.fwd.hh>
#include <protocols/hybridization/FoldTreeHybridize.hh>

#include <protocols/moves/Mover.hh>

#include <protocols/loops/Loops.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <utility/file/FileName.hh>
#include <utility/pointer/owning_ptr.hh>


#define cnt_hill_moves_without_increase 150
#define cnt_my_trial_frag_insertions 500

class CompleteAbinitioSampler : public protocols::abinitio::ClassicAbinitio {
public:
  CompleteAbinitioSampler(
			  core::fragment::FragSetCOP fragset_small,
			  core::fragment::FragSetCOP fragset_large,
			  core::kinematics::MoveMapCOP movemap
			  ) : ClassicAbinitio( fragset_small, fragset_large, movemap ) {};


  void apply(core::pose::Pose& pose_, FuncStats& stats);

  void apply( core::pose::Pose &pose ) override;
};


class StageRosettaSampler : public protocols::abinitio::ClassicAbinitio {
public:
  std::string rosetta_stage;
  protocols::simple_moves::ClassicFragmentMoverOP frag_mover_large;
  StageRosettaSampler(
		core::fragment::FragSetCOP fragset_small,
		core::fragment::FragSetCOP fragset_large,
		core::kinematics::MoveMapCOP movemap
		) : ClassicAbinitio( fragset_small, fragset_large, movemap ) {
    rosetta_stage = "stage4";
    frag_mover_large = protocols::simple_moves::ClassicFragmentMoverOP(  new protocols::simple_moves::ClassicFragmentMover(fragset_large, movemap) );
  };

  StageRosettaSampler(
		core::fragment::FragSetCOP fragset_small,
		core::fragment::FragSetCOP fragset_large,
		core::kinematics::MoveMapCOP movemap,
		std::string stage
		) : ClassicAbinitio( fragset_small, fragset_large, movemap ) {
    rosetta_stage = stage;
    frag_mover_large = protocols::simple_moves::ClassicFragmentMoverOP(  new protocols::simple_moves::ClassicFragmentMover(fragset_large, movemap) );
    set_skip_stage2(false);
  };

  void apply(core::pose::Pose& pose_, FuncStats& stats) ;
  void apply( core::pose::Pose &pose ) override;
};



class InitStagesSampler : public protocols::abinitio::ClassicAbinitio {
public:
  InitStagesSampler(
		    core::fragment::FragSetCOP fragset_small,
		    core::fragment::FragSetCOP fragset_large,
		    core::kinematics::MoveMapCOP movemap
		    ) : ClassicAbinitio( fragset_small, fragset_large, movemap ) {};

  void apply(core::pose::Pose& pose_, FuncStats& stats) ;
  void apply( core::pose::Pose &pose ) override;
};

class FragInsertionMover
{
public:
  protocols::simple_moves::ClassicFragmentMoverOP frag_mover, frag_mover_large;
  core::scoring::ScoreFunctionOP sfxn;

  FragInsertionMover(core::scoring::ScoreFunctionOP sfxn_, core::fragment::FragSetOP frag_set_, core::fragment::FragSetOP large_frag_set_) ;

  virtual void apply(core::pose::Pose& pose_, FuncStats& stats) ;

  virtual void apply(core::pose::Pose& pose_) ;
};





class PerturbFragMover : public FragInsertionMover {
public:
  PerturbFragMover(core::scoring::ScoreFunctionOP sfxn_, core::fragment::FragSetOP frag_set_, core::fragment::FragSetOP frag_set_large) : FragInsertionMover(sfxn_, frag_set_, frag_set_large) {
  }

  void apply(core::pose::Pose& pose_) ;
};


class LocalSearchFragMover : public FragInsertionMover {
public:
  LocalSearchFragMover(core::scoring::ScoreFunctionOP sfxn_, core::fragment::FragSetOP frag_set_, core::fragment::FragSetOP frag_set_large) : FragInsertionMover(sfxn_, frag_set_, frag_set_large) {
  }

  virtual void apply(core::pose::Pose& pose_, FuncStats& stats) ;

  void apply(core::pose::Pose& pose_) ;
};


class InitStagesMover : public FragInsertionMover {
public:
  core::scoring::ScoreFunctionOP sfxn;
  boost::shared_ptr<InitStagesSampler> sampler;
  core::fragment::FragSetCOP fragset_small, fragset_large;
  core::kinematics::MoveMapOP mm_;

  InitStagesMover(core::scoring::ScoreFunctionOP sfxn_, core::fragment::FragSetOP frag_set_, core::fragment::FragSetOP frag_set_large) : FragInsertionMover(sfxn_, frag_set_, frag_set_large) {
    mm_ = core::kinematics::MoveMapOP(new core::kinematics::MoveMap());
    mm_->set_bb(true);
    fragset_small = core::fragment::FragSetCOP(frag_set_);
    fragset_large = core::fragment::FragSetCOP(frag_set_large);
    sampler =  boost::shared_ptr<InitStagesSampler>(new InitStagesSampler(fragset_small, fragset_large, mm_->clone() ) );
  }

  void apply(core::pose::Pose& pose);

};

class CompleteAbinitioMover : public FragInsertionMover
{
public:
public:
  core::scoring::ScoreFunctionOP sfxn;
  boost::shared_ptr<CompleteAbinitioSampler> sampler;
  core::fragment::FragSetCOP fragset_small, fragset_large;
  core::kinematics::MoveMapOP mm_;

  CompleteAbinitioMover(core::scoring::ScoreFunctionOP sfxn_, core::fragment::FragSetOP frag_set_, core::fragment::FragSetOP frag_set_large) : FragInsertionMover(sfxn_, frag_set_, frag_set_large) {
    mm_ = core::kinematics::MoveMapOP(new core::kinematics::MoveMap());
    mm_->set_bb(true);
    fragset_small = core::fragment::FragSetCOP(frag_set_);
    fragset_large = core::fragment::FragSetCOP(frag_set_large);
    sampler =  boost::shared_ptr<CompleteAbinitioSampler>(new CompleteAbinitioSampler(fragset_small, fragset_large, mm_->clone() ) );
    sfxn = sfxn_;
  }

  void apply(core::pose::Pose& pose_);

};


class StageRosettaMover : public FragInsertionMover
{
public:
public:
  core::scoring::ScoreFunctionOP sfxn;
  boost::shared_ptr<StageRosettaSampler> sampler;
  core::fragment::FragSetCOP fragset_small, fragset_large;
  core::kinematics::MoveMapOP mm_;

  StageRosettaMover(core::scoring::ScoreFunctionOP sfxn_, core::fragment::FragSetOP frag_set_, core::fragment::FragSetOP frag_set_large, std::string stage = "stage4") : FragInsertionMover(sfxn_, frag_set_, frag_set_large) {
    mm_ = core::kinematics::MoveMapOP(new core::kinematics::MoveMap());
    mm_->set_bb(true);
    fragset_small = core::fragment::FragSetCOP(frag_set_);
    fragset_large = core::fragment::FragSetCOP(frag_set_large);
    sampler =  boost::shared_ptr<StageRosettaSampler>(new StageRosettaSampler(fragset_small, fragset_large, mm_->clone() , stage) );
    sfxn = sfxn_;
  }

  virtual void apply(core::pose::Pose& pose_);
};

class NoGreedyStageRosettaMover : public StageRosettaMover
{
public:
  NoGreedyStageRosettaMover(core::scoring::ScoreFunctionOP sfxn_, core::fragment::FragSetOP frag_set_, core::fragment::FragSetOP frag_set_large, std::string stage) : StageRosettaMover(sfxn_, frag_set_, frag_set_large, stage) {
  }

  void apply(core::pose::Pose& pose_, FuncStats& stats) ;
  void apply(core::pose::Pose& pose_) ;

};

class FragInsertionStrategy
{
public:
  enum FragMoverTypes {
    my_frag_insertion, greedy_search, no_greedy_search, stage_rosetta_mover, hybrid_mover, cartesian_mover
  };

  FragInsertionStrategy() {}

  struct FragOptions
  {
    core::scoring::ScoreFunctionOP scorefxn;
    core::fragment::FragSetOP frag_set_;
    core::fragment::FragSetOP frag_set_large;
    core::pose::PoseOP native_model;
    core::pose::PoseOP consen_model;
    std::string ss;
    std::string stage_name;
  };

  static boost::shared_ptr<FragInsertionMover> get(FragMoverTypes type, FragOptions opt ) {
    return get(type, opt.scorefxn, opt.frag_set_, opt.frag_set_large, opt.stage_name, opt.native_model);
  }

  static boost::shared_ptr<FragInsertionMover> get(FragMoverTypes type, core::scoring::ScoreFunctionOP scorefxn, core::fragment::FragSetOP frag_set_, core::fragment::FragSetOP frag_set_large, std::string stage_name, core::pose::PoseOP native_model) {

    // SUPER DANGER PLEASE FIX IT
     boost::shared_ptr<FragInsertionMover> frag_mover = boost::shared_ptr<FragInsertionMover>( new StageRosettaMover(scorefxn, frag_set_, frag_set_large, stage_name));
    return frag_mover;

    
    switch (type) {
    case my_frag_insertion: {
      frag_mover = boost::shared_ptr<FragInsertionMover>( new FragInsertionMover(scorefxn, frag_set_, frag_set_large));
      break;
    }
    case greedy_search: {
      frag_mover = boost::shared_ptr<FragInsertionMover>( new LocalSearchFragMover(scorefxn, frag_set_, frag_set_large));
      break;
    }
    case no_greedy_search: {
      frag_mover = boost::shared_ptr<FragInsertionMover>( new NoGreedyStageRosettaMover(scorefxn, frag_set_, frag_set_large, stage_name));
      break;
    }
    case stage_rosetta_mover: {
      frag_mover = boost::shared_ptr<FragInsertionMover>( new StageRosettaMover(scorefxn, frag_set_, frag_set_large, stage_name));
      break;
    }
    default:
      frag_mover = boost::shared_ptr<FragInsertionMover>( new FragInsertionMover(scorefxn, frag_set_, frag_set_large));
      break;
    }

    return frag_mover;
  }
};


typedef boost::shared_ptr<FragInsertionMover> FragInsertionMoverPtr;

#endif /* FRAGINSERTIONMOVER_H */
