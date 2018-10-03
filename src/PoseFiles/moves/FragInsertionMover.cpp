
#include "FitFunction.hpp"
#include "FragInsertionMover.hpp"

#include <core/scoring/Energies.hh>

void
CompleteAbinitioSampler::apply(core::pose::Pose& pose_, FuncStats& stats) {
  apply(pose_);
}

void
CompleteAbinitioSampler::apply( core::pose::Pose &pose ) {
  mc().reset_counters();

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
  if (rosetta_stage == "stage2") {
    stage2_cycles_ = 50;
    prepare_stage2( pose );
    do_stage2_cycles( pose );
    recover_low( pose , STAGE_2 );
    mc().reset_counters();
  }

  if (rosetta_stage == "stage3") {
    stage3_cycles_ = 15;
    prepare_stage3( pose );
    do_stage3_cycles( pose );
    recover_low( pose , STAGE_3b );
    mc().reset_counters();
  }

  if (rosetta_stage == "stage4") {
    stage4_cycles_ = 25;
    prepare_stage4( pose );
    do_stage4_cycles( pose );
    recover_low( pose , STAGE_4 );
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
  do
    {
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
    } while (cnt < 100);
}

void
FragInsertionMover::apply(core::pose::Pose& pose_) {
  core::pose::Pose inner_pose = pose_;
  double init_score = (*sfxn)(pose_);
  double current_score;
  int cnt = 0;
  do
    {
      frag_mover->apply(inner_pose);
      current_score = (*sfxn)(inner_pose);
      if (current_score < init_score) {
	pose_ = inner_pose;
	init_score = current_score;
      } else {
	inner_pose = pose_;
      }

      cnt++;
    } while (cnt < 50);
}


void
PerturbFragMover::apply(core::pose::Pose& pose_) {
  frag_mover->apply(pose_);
  (*sfxn)(pose_);
}

void
LocalSearchFragMover::apply(core::pose::Pose& pose_, FuncStats& stats) {
  core::pose::Pose inner_pose = pose_;
  double init_score = (*sfxn)(pose_);
  double current_score;
  int cnt = 0;

  int hill_moves_without_increase = 0;
  do
    {
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

    } while (hill_moves_without_increase < 50);
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
    }while (hill_moves_without_increase < 50);
  //} while (cnt < 500);
}


void
ILSFragMover::apply(core::pose::Pose& pose_) {
  double init_score = (*sfxn)(pose_);
  core::pose::Pose backup_pose_ = pose_;
  perturb_fragmover->apply(pose_);
  ls_fragmover->apply(pose_);

  if (init_score < (*sfxn)(pose_)) {
    pose_ = backup_pose_;
  }

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
NoGreedyStageRosettaMover::apply(core::pose::Pose& pose_) {
  core::pose::Pose pose_backup = pose_;
  (*sfxn)(pose_);
  sampler->init(pose_);
  sampler->apply(pose_);
  (*sfxn)(pose_);
}
