
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


class HybridRosettaInserter : public FragInsertionMover {
public:

  class HybridRosettaSampler : public protocols::hybridization::HybridizeProtocol {
  public:
    HybridRosettaSampler(std::string stage,
			 core::fragment::FragSetCOP fragset_small,
			 core::fragment::FragSetCOP fragset_large,
			 core::kinematics::MoveMapCOP movemap,
			 utility::vector1 <core::pose::PoseOP> templates_in,
			 utility::vector1 <core::Real> template_weights_in,
			 core::scoring::ScoreFunctionOP stage1_scorefxn_in,
			 core::scoring::ScoreFunctionOP stage2_scorefxn_in,
			 core::scoring::ScoreFunctionOP fa_scorefxn_in,
			 std::string frag3_fn,
			 std::string frag9_fn,
			 std::string cen_cst_in,
			 std::string fa_cst_in
			 ) : protocols::hybridization::HybridizeProtocol(
									 templates_in,
									 template_weights_in,
									 stage1_scorefxn_in,
									 stage2_scorefxn_in,
									 fa_scorefxn_in,
									 frag3_fn,
									 frag9_fn,
									 cen_cst_in,
									 fa_cst_in
									 ) {
      if (stage =="stage1") {
	setup_cycles(stage2);
      }
      if (stage =="stage2") {
	setup_cycles(stage3);
      }
      if (stage =="stage3") {
	setup_cycles(stage3);
      }
      if (stage =="stage4") {
	setup_cycles(stage4);
      }
      if (stage =="complete") {
	setup_cycles(complete);
      }
    }

    enum stages_type {
      stage1, stage2, stage3, stage4, complete
    };

    stages_type current_stage;

    void apply( core::pose::Pose & pose );
  public:
    void setup_cycles(stages_type stage_num) {
      stage1_increase_cycles_ = 1.0;
      min_after_stage1_ = false;
      if (stage_num == stage1) {
	stage1_1_cycles_    = 1;
	stage1_2_cycles_    = 0;
	stage1_3_cycles_    = 0;
	stage1_4_cycles_    = 0;
      }

      if (stage_num == stage2) {
	stage1_1_cycles_    = 0;
	stage1_2_cycles_    = 1;
	stage1_3_cycles_    = 0;
	stage1_4_cycles_    = 0;

      }
      if (stage_num == stage3) {
	stage1_1_cycles_    = 0;
	stage1_2_cycles_    = 0;
	stage1_3_cycles_    = 1;
	stage1_4_cycles_    = 0;

      }
      if (stage_num == stage4) {
	stage1_1_cycles_    = 0;
	stage1_2_cycles_    = 0;
	stage1_3_cycles_    = 0;
	stage1_4_cycles_    = 1;
      }

      if (stage_num == complete) {
	min_after_stage1_ = true;
	//stage1_1_cycles_    = 0;
	//stage1_2_cycles_    = 0;
	//stage1_3_cycles_    = 0;
	//stage1_4_cycles_    = 30;
      }
    }
  };

  // members
  boost::shared_ptr<HybridRosettaSampler> hybrid_protocol;
  core::kinematics::MoveMapOP mm_;
  std::string rosetta_stage;

  // constructor
  HybridRosettaInserter(core::scoring::ScoreFunctionOP sfxn_, core::fragment::FragSetOP frag_set_, core::fragment::FragSetOP large_frag_set_) : FragInsertionMover(sfxn_, frag_set_, large_frag_set_) {
  }


  // builder
  static
  boost::shared_ptr<FragInsertionMover> Builder(std::string stage, core::scoring::ScoreFunctionOP sfxn_,  core::fragment::FragSetOP fragset_small, core::fragment::FragSetOP fragset_large, core::pose::PoseOP model ) {
    utility::vector1 <core::pose::PoseOP> templates_in;
    templates_in.push_back(model);
    utility::vector1 <core::Real> template_weights_in;
    template_weights_in.push_back(1.0);
    core::scoring::ScoreFunctionOP stage1_scorefxn_in = sfxn_;
    //core::scoring::ScoreFunctionOP stage2_scorefxn_in = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score4_smooth_cart").c_str());
    core::scoring::ScoreFunctionOP stage2_scorefxn_in = sfxn_;
    core::scoring::ScoreFunctionOP fa_scorefxn_in = core::scoring::ScoreFunctionFactory::create_score_function(std::string("talaris2014_cart").c_str());
    std::string frag3_fn = "./input_files/info_1wit/boinc_vf_aa1wit_03_05.200_v1_3";
    std::string frag9_fn = "./input_files/info_1wit/boinc_vf_aa1wit_09_05.200_v1_3";
    std::string cen_cst_in = "AUTO";
    std::string fa_cst_in = "AUTO";
    core::kinematics::MoveMapOP movemap;
    movemap = core::kinematics::MoveMapOP(new core::kinematics::MoveMap());
    movemap->set_bb(true);
    boost::shared_ptr<HybridRosettaSampler> protocol_ = boost::shared_ptr<HybridRosettaSampler>(new HybridRosettaSampler( stage,
fragset_small, fragset_large, movemap,	templates_in, template_weights_in,
stage1_scorefxn_in,stage2_scorefxn_in,fa_scorefxn_in,frag3_fn,frag9_fn,cen_cst_in,fa_cst_in
));
    protocol_->set_stage1_increase_cycles(1.0);
    protocol_->set_stage2_increase_cycles(0);
    protocol_->set_batch_relax(0);
    boost::shared_ptr<HybridRosettaInserter> frag_inserter = boost::shared_ptr<HybridRosettaInserter>(new HybridRosettaInserter(sfxn_, fragset_small, fragset_large));
    frag_inserter->rosetta_stage = stage;
    frag_inserter->hybrid_protocol = protocol_;
    return frag_inserter;
  }

  void apply(core::pose::Pose& pose_, FuncStats& stats) {
    apply(pose_);
  }

  void apply( core::pose::Pose &pose ) override;

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

  void apply(core::pose::Pose& pose_, FuncStats& stats) ;
  void apply(core::pose::Pose& pose_) ;

};

class FragInsertionStrategy
{
public:
  enum FragMoverTypes {
    my_frag_insertion, greedy_search, no_greedy_search, stage_rosetta_mover, hybrid_mover, ILS_as_julia
  };

  FragInsertionStrategy() {}

  struct FragOptions
  {
    core::scoring::ScoreFunctionOP scorefxn;
    core::fragment::FragSetOP frag_set_;
    core::fragment::FragSetOP frag_set_large;
    core::pose::PoseOP native_model;
    std::string stage_name;
  };

  static boost::shared_ptr<FragInsertionMover> get(FragMoverTypes type, FragOptions opt ) {
    return get(type, opt.scorefxn, opt.frag_set_, opt.frag_set_large, opt.stage_name, opt.native_model);
  }

  static boost::shared_ptr<FragInsertionMover> get(FragMoverTypes type, core::scoring::ScoreFunctionOP scorefxn, core::fragment::FragSetOP frag_set_, core::fragment::FragSetOP frag_set_large, std::string stage_name, core::pose::PoseOP native_model) {
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
    case no_greedy_search: {
      frag_mover = boost::shared_ptr<FragInsertionMover>( new NoGreedyStageRosettaMover(scorefxn, frag_set_, frag_set_large, stage_name));
      break;
    }
    case hybrid_mover: {
      frag_mover = HybridRosettaInserter::Builder(stage_name, scorefxn, frag_set_, frag_set_large, native_model);
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
