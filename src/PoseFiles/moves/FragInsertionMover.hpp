
#ifndef FRAGINSERTIONMOVER_H
#define FRAGINSERTIONMOVER_H

#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MonteCarlo.hh>
#include "FitFunction.hpp"

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
  StageRosettaSampler(
		core::fragment::FragSetCOP fragset_small,
		core::fragment::FragSetCOP fragset_large,
		core::kinematics::MoveMapCOP movemap
		) : ClassicAbinitio( fragset_small, fragset_large, movemap ) {
    rosetta_stage = "stage4";
  };

  StageRosettaSampler(
		core::fragment::FragSetCOP fragset_small,
		core::fragment::FragSetCOP fragset_large,
		core::kinematics::MoveMapCOP movemap,
		std::string stage
		) : ClassicAbinitio( fragset_small, fragset_large, movemap ) {
    rosetta_stage = stage;
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

class ILSFragMover : public FragInsertionMover {
public:
  boost::shared_ptr<LocalSearchFragMover> ls_fragmover;
  boost::shared_ptr<PerturbFragMover> perturb_fragmover;

  ILSFragMover(core::scoring::ScoreFunctionOP sfxn_, core::fragment::FragSetOP frag_set_, core::fragment::FragSetOP frag_set_large) : FragInsertionMover(sfxn_, frag_set_, frag_set_large) {
    ls_fragmover = boost::shared_ptr<LocalSearchFragMover>(new LocalSearchFragMover(sfxn_, frag_set_, frag_set_large));
    perturb_fragmover = boost::shared_ptr<PerturbFragMover>(new PerturbFragMover(sfxn_, frag_set_, frag_set_large));
  }


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

  void apply(core::pose::Pose& pose_) ;

};

class FragInsertionStrategy
{
public:
  enum FragMoverTypes {
    my_frag_insertion, greedy_search, stage_rosetta_mover, ILS_as_julia
  };

  FragInsertionStrategy() {}

  struct FragOptions
  {
    core::scoring::ScoreFunctionOP scorefxn;
    core::fragment::FragSetOP frag_set_;
    core::fragment::FragSetOP frag_set_large;
    std::string stage_name;
  };

  static boost::shared_ptr<FragInsertionMover> get(FragMoverTypes type, FragOptions opt ) {
    return get(type, opt.scorefxn, opt.frag_set_, opt.frag_set_large, opt.stage_name);
  }

  static boost::shared_ptr<FragInsertionMover> get(FragMoverTypes type, core::scoring::ScoreFunctionOP scorefxn, core::fragment::FragSetOP frag_set_, core::fragment::FragSetOP frag_set_large, std::string stage_name) {

    boost::shared_ptr<FragInsertionMover> frag_mover = boost::shared_ptr<FragInsertionMover>( new ILSFragMover(scorefxn, frag_set_, frag_set_large));
    switch (type) {
    case my_frag_insertion: {
      frag_mover = boost::shared_ptr<FragInsertionMover>( new FragInsertionMover(scorefxn, frag_set_, frag_set_large));
      break;
    }
    case greedy_search: {
      frag_mover = boost::shared_ptr<FragInsertionMover>( new LocalSearchFragMover(scorefxn, frag_set_, frag_set_large));
      break;
    }
    case stage_rosetta_mover: {
      frag_mover = boost::shared_ptr<FragInsertionMover>( new StageRosettaMover(scorefxn, frag_set_, frag_set_large, stage_name));
      break;
    }
    case ILS_as_julia: {
      frag_mover = boost::shared_ptr<FragInsertionMover>( new ILSFragMover(scorefxn, frag_set_, frag_set_large));
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
