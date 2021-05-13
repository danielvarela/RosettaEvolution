
#include "FitFunction.hpp"
#include "FragInsertionMover.hpp"

#include <protocols/hybridization/FoldTreeHybridize.fwd.hh>
#include <core/scoring/Energies.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>

#include <core/init/score_function_corrections.hh>

#include <protocols/hybridization/HybridizeProtocolCreator.hh>
#include <protocols/hybridization/HybridizeProtocol.hh>
#include <protocols/hybridization/FoldTreeHybridize.hh>
#include <protocols/hybridization/CartesianHybridize.hh>
#include <protocols/hybridization/TemplateHistory.fwd.hh>
#include <protocols/hybridization/util.hh>
#include <protocols/hybridization/DomainAssembly.hh>
#include <protocols/hybridization/DDomainParse.hh>
#include <protocols/hybridization/TMalign.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/minimization_packing/MinMover.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/util.hh>

#include <protocols/loops/util.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>

#include <protocols/rigid/RB_geometry.hh>

#include <protocols/electron_density/SetupForDensityScoringMover.hh>

// dynamic fragpick
#include <protocols/moves/DsspMover.hh>
#include <core/fragment/picking_old/vall/util.hh>

#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentStruct.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/CrystInfo.hh>
#include <core/import_pose/import_pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragID_Iterator.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/util.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/SWAligner.hh>
#include <core/sequence/SimpleScoringScheme.hh>
#include <core/sequence/ScoringScheme.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/func/CircularSplineFunc.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// task operation
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

#include <basic/datacache/DataMap.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// utility
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <numeric/random/WeightedSampler.hh>
#include <ObjexxFCL/format.hh>
#include <boost/foreach.hpp>

#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>

// evaluation
#include <core/scoring/rms_util.hh>
#include <protocols/comparative_modeling/coord_util.hh>

// option
#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh> // strand pairings
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh> // check pre talaris

//docking
#include <protocols/docking/DockingLowRes.hh>
#include <protocols/symmetric_docking/SymDockingLowRes.hh>
#include <protocols/minimization_packing/symmetry/SymMinMover.hh>

#include <string>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


static basic::Tracer TR( "protocols.hybridization.HybridizeProtocol" );


void
CompleteAbinitioSampler::apply(core::pose::Pose& pose_, FuncStats& stats) {
  apply(pose_);
}

void
CompleteAbinitioSampler::apply( core::pose::Pose &pose ) {
  mc().reset_counters();
  stage1_cycles_ = stage1_cycles_ * 10.0;
  stage2_cycles_ = stage2_cycles_ * 10.0;
  stage3_cycles_ = stage3_cycles_ * 10.0;
  stage4_cycles_ = stage4_cycles_ * 10.0;
  prepare_stage1( pose );
  do_stage1_cycles( pose );
  recover_low( pose , STAGE_1 );
  // mc().show_counters();
  mc().reset_counters();
  //current_scorefxn().show( std::cout, pose);

  prepare_stage2( pose );
  do_stage2_cycles( pose );
  recover_low( pose , STAGE_2 );
  //mc().show_counters();
  mc().reset_counters();
  //current_scorefxn().show( std::cout, pose);

  prepare_stage3( pose );
  do_stage3_cycles( pose );
  recover_low( pose , STAGE_3b );
  // mc().show_counters();
  mc().reset_counters();
  //current_scorefxn().show( std::cout, pose);

  prepare_stage4( pose );
  do_stage4_cycles( pose );
  recover_low( pose , STAGE_4 );
  //    mc().show_counters();
  mc().reset_counters();
  // current_scorefxn().show( std::cout, pose);
}


void
StageRosettaSampler::apply(core::pose::Pose& pose_, FuncStats& stats) {
  apply(pose_);
}

void
StageRosettaSampler::apply( core::pose::Pose &pose ) {
  using namespace protocols::moves;
  /* stage2_cycles_ = 25;
  stage3_cycles_ = 15;
  stage4_cycles_ = 25;
 */ 

  stage2_cycles_ = stage2_cycles_ * 0.01;
  stage3_cycles_ = stage3_cycles_ * 0.01;
  stage4_cycles_ = stage4_cycles_ * 0.01;

  if (rosetta_stage == "stage2") {
    prepare_stage2( pose );
    do_stage2_cycles( pose );
    recover_low( pose , STAGE_2 );
    mc().show_state();
    mc().show_counters();
    mc().reset_counters();
  }

  if (rosetta_stage == "stage3") {
    prepare_stage3( pose );
    do_stage3_cycles( pose );
    mc().reset_counters();
  }

  if (rosetta_stage == "stage4") {
    prepare_stage4( pose );
    //mc().set_temperature(1.0);
    do_stage4_cycles( pose );
    //recover_low( pose , STAGE_4 );
    mc().show_counters();
    mc().reset_counters();
  }
}


void
InitStagesSampler::apply(core::pose::Pose& pose_, FuncStats& stats) {
  apply(pose_);
}

void
InitStagesSampler::apply( core::pose::Pose &pose ) {
  //  stage2_cycles_ = 100;
  bSkipStage2_ = true;
  bSkipStage3_ = true;
  bSkipStage4_ = true;

  // mc().set_autotemp(true, 2.0);
  // mc().set_temperature(2.0);
  // mc().reset(pose);
  prepare_stage1( pose );
  do_stage1_cycles( pose );
  recover_low( pose, STAGE_1 );
  // prepare_stage2( pose );
  // do_stage2_cycles( pose );

}


FragInsertionMover::FragInsertionMover(core::scoring::ScoreFunctionOP sfxn_, core::fragment::FragSetOP frag_set_, core::fragment::FragSetOP large_frag_set_) {
  core::kinematics::MoveMapOP mm_ = core::kinematics::MoveMapOP(new core::kinematics::MoveMap());
  mm_->set_bb(true);

  //    std::cout << "frag_mover op" << std::endl;
  frag_mover = protocols::simple_moves::ClassicFragmentMoverOP( new protocols::simple_moves::ClassicFragmentMover(frag_set_, mm_) );
  // std::cout << "frag_mover op 2" << std::endl;
  frag_mover_large = protocols::simple_moves::ClassicFragmentMoverOP(  new protocols::simple_moves::ClassicFragmentMover(large_frag_set_, mm_) );
  //std::cout << "finish " << std::endl;
  sfxn = sfxn_;
}

void
FragInsertionMover::apply(core::pose::Pose& pose_, FuncStats& stats) {
  core::pose::Pose inner_pose = pose_;
  double init_score = (*sfxn)(pose_);
  double current_score;
  int cnt = 0;
  do {
      frag_mover->apply(inner_pose);
      current_score = (*sfxn)(inner_pose);
      if (stats.find("total_tries") != stats.end()) {
	stats["total_tries"]++;
      } else {
	stats["total_tries"] = 1;
      }

      if (current_score < init_score) {
	if (stats.find("accepted_fragments") != stats.end()) {
	  stats["accepted_fragments"]++;
	} else {
	  stats["accepted_fragments"] = 1;
	}

	pose_ = inner_pose;
	init_score = current_score;
      } else {
	inner_pose = pose_;
      }

      cnt++;
    } while (cnt < cnt_my_trial_frag_insertions);
}

void
FragInsertionMover::apply(core::pose::Pose& pose_) {
  core::pose::Pose inner_pose = pose_;
  double init_score = (*sfxn)(pose_);
  double current_score;
  int cnt = 0;
  do  {
      frag_mover->apply(inner_pose);
      current_score = (*sfxn)(inner_pose);
      if (current_score < init_score) {
	pose_ = inner_pose;
	init_score = current_score;
      } else {
	inner_pose = pose_;
      }

      cnt++;
    } while (cnt < cnt_my_trial_frag_insertions);
}


void
PerturbFragMover::apply(core::pose::Pose& pose_) {
  frag_mover_large->apply(pose_);
  //  frag_mover->apply(pose_);
  double score = (*sfxn)(pose_);
  std::cout << "pose_ " <<  score << std::endl;
}

void
LocalSearchFragMover::apply(core::pose::Pose& pose_, FuncStats& stats) {
  core::pose::Pose inner_pose = pose_;
  double init_score = (*sfxn)(pose_);
  double current_score;
  int cnt = 0;

  int hill_moves_without_increase = 0;
  do {
      frag_mover->apply(inner_pose);
      current_score = (*sfxn)(inner_pose);
      if (stats.find("total_tries") != stats.end()) {
	stats["total_tries"]++;
      } else {
	stats["total_tries"] = 1;
      }


      if (current_score < init_score) {
	if (stats.find("accepted_fragments") != stats.end()) {
	  stats["accepted_fragments"]++;
	} else {
	  stats["accepted_fragments"] = 1;
	}

	pose_ = inner_pose;
	init_score = current_score;
      } else {
	inner_pose = pose_;
	hill_moves_without_increase++;
      }

      cnt++;
      // } while (cnt < 50);

    } while (hill_moves_without_increase < 150);
}

void
LocalSearchFragMover::apply(core::pose::Pose& pose_) {
  core::pose::Pose inner_pose = pose_;
  double init_score = (*sfxn)(pose_);
  double current_score;
  int cnt = 0;

  int hill_moves_without_increase = 0;
  do {
      frag_mover->apply(inner_pose);
      current_score = (*sfxn)(inner_pose);

      if (current_score < init_score) {
	pose_ = inner_pose;
	init_score = current_score;
      } else {
	inner_pose = pose_;
	hill_moves_without_increase++;
      }

      cnt++;
    }while (hill_moves_without_increase < 150);
  //} while (cnt < 500);
}




void
InitStagesMover::apply(core::pose::Pose& pose) {
  sampler->init(pose);
  sampler->apply(pose);
}


void
CompleteAbinitioMover::apply(core::pose::Pose& pose_) {
  core::pose::Pose pose_backup = pose_;
  double init_score = (*sfxn)(pose_);
  sampler->init(pose_);
  sampler->apply(pose_);
  std::cout << "ENTRA EN COMPLETE AB INITIO" << std::endl;
  exit(1);
  double final_score = (*sfxn)(pose_);
  if (init_score < final_score) {
    pose_ = pose_backup;
    (*sfxn)(pose_);
  }
}

void
StageRosettaMover::apply(core::pose::Pose& pose_) {
  core::pose::Pose pose_backup = pose_;
  double init_score = (*sfxn)(pose_);
  sampler->init(pose_);
  sampler->apply(pose_);

  double final_score = (*sfxn)(pose_);
  if (init_score < final_score) {
    pose_ = pose_backup;
    (*sfxn)(pose_);
  }
}

void
NoGreedyStageRosettaMover::apply(core::pose::Pose& pose_, FuncStats& stats) {
  apply(pose_);
}

void
NoGreedyStageRosettaMover::apply(core::pose::Pose& pose_) {
  //core::pose::Pose pose_backup = pose_;
  double init_score = (*sfxn)(pose_);
  sampler->init(pose_);
  sampler->apply(pose_);
  double score = (*sfxn)(pose_);
}

