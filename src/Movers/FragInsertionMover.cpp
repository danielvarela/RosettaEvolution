
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
#include <protocols/symmetry/SetupForSymmetryMover.hh>
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

  stage2_cycles_ = stage2_cycles_ * 0.1;
  stage3_cycles_ = stage3_cycles_ * 0.1;
  stage4_cycles_ = stage4_cycles_ * 0.1;

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

void
HybridRosettaInserter::apply(core::pose::Pose &pose ) {
  //std::cout << "Hybrid apply" << std::endl;
  hybrid_protocol->apply(pose);
}


void
CartesianRosettaInserter::apply(core::pose::Pose &pose ) {
  //std::cout << "Hybrid apply" << std::endl;
  cartesian_protocol->apply(pose);
  double score = (*sfxn)(pose);
}

void HybridRosettaInserter::HybridRosettaSampler::apply( core::pose::Pose & pose )
{
  using namespace protocols::hybridization; //::HybridizeProtocol::apply(pose);
  // using TemplateHistoryOP = typename protocols::hybridization::TemplateHistory::TemplateHistoryOP;
  // using FoldTreeHybridize = typename protocols::hybridization::FoldTreeHybridize;
  // typedef boost::shared_ptr<FoldTreeHybridize> FoldTreeHybridizeOP;
  // using FoldTreeHybridize = typename protocols::hybridization::FoldTreeHybridize;
  // using FoldTreeHybridizeOP = typename protocols::hybridization::FoldTreeHybridizeOP;
  using namespace core;
  using namespace kinematics;
  using namespace sequence;
  using namespace pack;
  using namespace task;
  using namespace operation;
  using namespace scoring;
  using namespace constraints;
  using namespace protocols::loops;
  using namespace protocols::moves;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace core::pose::datacache;
  using namespace core::io::silent;
  using namespace ObjexxFCL::format;

  // number of residues in asu without VRTs
  core::Size nres_tgt = pose.size();
  core::Size nres_protein_tgt = pose.size();
  core::conformation::symmetry::SymmetryInfoCOP symm_info;
  if ( core::pose::symmetry::is_symmetric(pose) ) {
    auto & SymmConf (
		     dynamic_cast<core::conformation::symmetry::SymmetricConformation &> ( pose.conformation()) );
    symm_info = SymmConf.Symmetry_Info();
    nres_tgt = symm_info->num_independent_residues();
    nres_protein_tgt = symm_info->num_independent_residues();
  }

  if ( pose.residue(nres_tgt).aa() == core::chemical::aa_vrt ) nres_tgt--;
  while ( !pose.residue(nres_protein_tgt).is_protein() ) nres_protein_tgt--;

  //save necessary constraint in pose
  core::scoring::constraints::ConstraintSetOP save_pose_constraint_set( new core::scoring::constraints::ConstraintSet() ) ;
  if ( keep_pose_constraint_ ) {
    save_pose_constraint_set = pose.constraint_set()->clone();
    core::scoring::constraints::remove_nonbb_constraints(pose);
  }

  if ( pose.is_fullatom() ) {
    protocols::moves::MoverOP tocen( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );
    tocen->apply(pose);
  }

  // make fragments if we don't have them at this point
  check_and_create_fragments( pose );

  // domain parse templates if we do not have domain definitions
  domain_parse_templates(nres_tgt);

  // select random fragment sets if >2 are provided
  core::fragment::FragSetOP frags_small = fragments_small_[numeric::random::rg().random_range(1,fragments_small_.size())];
  core::fragment::FragSetOP frags_big = fragments_big_[numeric::random::rg().random_range(1,fragments_big_.size())];
  TR.Info << "FRAGMENTS small max length: " << frags_small->max_frag_length() << std::endl;
  TR.Info << "FRAGMENTS big max length: " << frags_big->max_frag_length() << std::endl;

  // starting structure pool (for batch relax)
  std::vector < SilentStructOP > post_centroid_structs;
  bool need_more_samples = true;

  // main centroid generation loop
  core::Real gdtmm = 0.0;
  //while ( need_more_samples ) {
    need_more_samples = false;
    // pick starting template
    core::Size initial_template_index = pick_starting_template();
    TR << "Using initial template: " << I(4,initial_template_index) << " " << template_fn_[initial_template_index] << std::endl;

    // ensure
    //    1)that no CONTIGS cross multiple CHAINS
    //    2)that every input CHAIN has a chunk in it
    // if not, add the corresponding contig as a chunk
    utility::vector1< int > cuts = pose.fold_tree().cutpoints();
    protocols::loops::Loops const &contigs_old = template_contigs_[initial_template_index];
    protocols::loops::Loops contigs_new;
    for ( core::Size j=1; j<=contigs_old.num_loop(); ++j ) {
      protocols::loops::Loop contig_j = contigs_old[j];

      for ( core::Size i=1; i<=cuts.size(); ++i ) {
	int nextcut = cuts[i];
	if ( nextcut > (int)nres_protein_tgt ) break;

	core::Size start = templates_[initial_template_index]->pdb_info()->number(contig_j.start());
	core::Size stop  = templates_[initial_template_index]->pdb_info()->number(contig_j.stop());

	if ( (int)start <= nextcut && (int)stop > nextcut ) {
	  // find template res corresponding to cut
	  for ( core::Size k=contig_j.start(); k<=contig_j.stop(); ++k ) {
	    if ( templates_[initial_template_index]->pdb_info()->number(k) == nextcut ) {
	      // ensure contig >= 3 res
	      if ( k > contig_j.start()+1 ) {
		contigs_new.push_back( protocols::loops::Loop(contig_j.start(), k) );
	      }
	      TR << "Split contig ("<<contig_j.start()<<","<<contig_j.stop()<< ") at " << k << " [pdbnum: "<<start<<"/"<<stop<<"/" << nextcut << "]" << std::endl;
	      contig_j.set_start( k+1 );
	    }
	  }

	}
      }
      // ensure contig >= 3 res
      if ( contig_j.stop() > contig_j.start()+1 ) {
	contigs_new.push_back( contig_j );
      }
    }
    template_contigs_[initial_template_index] = contigs_new;

    protocols::loops::Loops &chunks = template_chunks_[initial_template_index];
    for ( core::Size i=0; i<=cuts.size(); ++i ) {
      int prevcut = (i==0) ? 1 : cuts[i];
      int nextcut = (i==cuts.size()) ? nres_protein_tgt : cuts[i+1];
      if ( nextcut > (int)nres_protein_tgt ) break;
      bool haschunk = false;
      for ( core::Size j=1; j<=chunks.num_loop() && !haschunk; ++j ) {
	protocols::loops::Loop const &chunk_j = chunks[j];
	core::Size start = templates_[initial_template_index]->pdb_info()->number(chunk_j.start());
	core::Size stop  = templates_[initial_template_index]->pdb_info()->number(chunk_j.stop());
	if ( (int)start > prevcut && (int)stop <= nextcut ) haschunk=true;
      }

      if ( ! haschunk ) {
	TR << "Segment ("<<prevcut<<","<<nextcut<<") has no chunks!  Adding contigs!" << std::endl;
	bool hascontig=false;
	for ( core::Size j=1; j<=contigs_new.num_loop(); ++j ) {
	  protocols::loops::Loop const &contig_j = contigs_new[j];
	  core::Size start = templates_[initial_template_index]->pdb_info()->number(contig_j.start());
	  core::Size stop  = templates_[initial_template_index]->pdb_info()->number(contig_j.stop());
	  if ( (int)start > prevcut && (int)stop <= nextcut ) {
	    hascontig = true;
	    chunks.push_back( contig_j );
	  }
	}
	if ( !hascontig ) {
	  TR << "Warning!  No contigs found for segment!" << std::endl;
	  TR << "If you are not using the 'randomize=X' option there is likely something wrong with your input and models will not be reasonable!" << std::endl;
	}
      }
    }

    // (0) randomize chains
    if ( randomize_chains_[initial_template_index].size() > 0 ) {
      runtime_assert ( templates_[initial_template_index]->pdb_info() );
      numeric::xyzVector< core::Real > comFixed(0,0,0), comMoving(0,0,0);
      core::Real nFixed=0, nMoving=0;

      TR << "Randomize >" << randomize_chains_[initial_template_index].size() << "<" << std::endl;
      for ( core::Size q=1; q<=randomize_chains_[initial_template_index].size(); ++q ) {
	TR << "Randomize >" << randomize_chains_[initial_template_index][q] << "<" << std::endl;
      }

      for ( core::Size i=1; i<=templates_[initial_template_index]->size(); ++i ) {
	char chain = templates_[initial_template_index]->pdb_info()->chain(i);
	if ( std::find(
		       randomize_chains_[initial_template_index].begin(),
		       randomize_chains_[initial_template_index].end(),
		       chain) != randomize_chains_[initial_template_index].end()
	     ) {
	  comMoving += templates_[initial_template_index]->residue(i).xyz(2);
	  nMoving++;
	} else {
	  comFixed += templates_[initial_template_index]->residue(i).xyz(2);
	  nFixed++;
	}
      }

      // sanity check
      if ( nMoving == 0 ) {
	utility_exit_with_message("Randomize chain enabled but chain(s) not found!");
      }

      comMoving /= nMoving;
      comFixed /= nFixed;

      // randomize orientation
      numeric::xyzMatrix< core::Real > R = protocols::geometry::random_reorientation_matrix( 360, 360 );

      // randomize offset
      numeric::xyzVector< core::Real > T = numeric::random::random_point_on_unit_sphere< core::Real >( numeric::random::rg() );

      // perturb & slide into contact
      bool done=false, compacting=false;
      core::Real DIST = 100.0;
      core::scoring::ScoreFunction sf_vdw;
      core::Real vdw_score = -1.0, step_size = -2.0, min_step_size = 0.25;
      sf_vdw.set_weight( core::scoring::vdw, 1.0 );
      core::pose::Pose template_orig = *(templates_[initial_template_index]);
      while ( !done && DIST > 0 ) {
	utility::vector1< id::AtomID > atm_ids;
	utility::vector1< numeric::xyzVector< core::Real> > atm_xyzs;

	for ( core::Size i=1; i<=templates_[initial_template_index]->size(); ++i ) {
	  core::conformation::Residue const & rsd_i = template_orig.residue(i);
	  char chain = templates_[initial_template_index]->pdb_info()->chain(i);
	  if ( std::find(
			 randomize_chains_[initial_template_index].begin(),
			 randomize_chains_[initial_template_index].end(),
			 chain) != randomize_chains_[initial_template_index].end()
	       ) {
	    for ( core::Size j = 1; j <= rsd_i.natoms(); ++j ) {
	      id::AtomID id( j,i );
	      atm_ids.push_back( id );
	      numeric::xyzVector< core::Real > new_ij = R*(rsd_i.xyz(j) - comMoving) + comFixed + DIST * T;
	      atm_xyzs.push_back( new_ij );
	    }
	  }
	}
	templates_[initial_template_index]->batch_set_xyz( atm_ids, atm_xyzs );
	core::Real score_d = sf_vdw( *(templates_[initial_template_index]) );
	if ( vdw_score < 0 ) vdw_score = score_d;

	if ( !compacting && score_d > vdw_score + 2.0 ) {
	  compacting = true;
	  step_size = min_step_size;
	} else if ( compacting  && vdw_score <= score_d + 2.0 ) {
	  done = true;
	}

	TR << "D=" << DIST << " score=" << score_d << std::endl;

	DIST += step_size;
      }
    }

    // (1) steal hetatms from template
    // utility::vector1< std::pair< core::Size,core::Size > > hetatms;
    // if ( add_hetatm_ ) {
    //   for ( Size ires=1; ires <= templates_[initial_template_index]->size(); ++ires ) {
    // 	if ( templates_[initial_template_index]->pdb_info()->number(ires) > (int)nres_tgt ) {
    // 	  TR.Debug << "Insert hetero residue: " << templates_[initial_template_index]->residue(ires).name3() << std::endl;
    // 	  if ( templates_[initial_template_index]->residue(ires).is_polymer()
    // 	       && !templates_[initial_template_index]->residue(ires).is_lower_terminus()
    // 	       && !pose.residue(pose.size()).is_upper_terminus() ) {
    // 	    pose.append_residue_by_bond(templates_[initial_template_index]->residue(ires));
    // 	  } else {
    // 	    pose.append_residue_by_jump(templates_[initial_template_index]->residue(ires), 1);
    // 	  }
    // 	  hetatms.push_back( std::make_pair( ires, pose.size() ) );
    // 	}
    //   }
    //   pose.conformation().chains_from_termini();
    // }

    // (2) realign structures per-domain
    if ( realign_domains_ ) {
      //
      if ( std::find( non_null_template_indices_.begin(), non_null_template_indices_.end(), initial_template_index )
	   != non_null_template_indices_.end() ) {
	align_templates_by_domain(templates_[initial_template_index], domains_all_templ_[initial_template_index]);
      } else {
	if ( non_null_template_indices_.size() == 0 ) {
	  TR << "No non-extended templates.  Skipping alignment." << std::endl;
	} else {
	  TR << "Cannot align to extended template! Choosing a random template instead." << std::endl;
	  core::Size aln_target = numeric::random::random_range(1, non_null_template_indices_.size() );
	  align_templates_by_domain(templates_[aln_target], domains_all_templ_[aln_target]);
	}
      }
    }

    // (3) apply symmetry
    std::string symmdef_file = symmdef_files_[initial_template_index];
    if ( !symmdef_file.empty() && symmdef_file != "NULL" ) {
      protocols::symmetry::SetupForSymmetryMover makeSymm( symmdef_file );
      makeSymm.apply(pose);

      //fpd   to get the right rotamer set we need to do this
      basic::options::option[basic::options::OptionKeys::symmetry::symmetry_definition].value( "dummy" );

      // xyz copy hetatms (properly handle cases where scoring subunit is not the first)
      // if ( add_hetatm_ ) {
      // 	for ( Size ihet=1; ihet <= hetatms.size(); ++ihet ) {
      // 	  core::conformation::Residue const &res_in = templates_[initial_template_index]->residue(hetatms[ihet].first);
      // 	  for ( Size iatm=1; iatm<=res_in.natoms(); ++iatm ) {
      // 	    core::id::AtomID tgt(iatm,hetatms[ihet].second);
      // 	    pose.set_xyz( tgt, res_in.xyz( iatm ) );
      // 	  }
      // 	}
      // }

    }

    // (4) add a virtual if we have a map or coord csts
    //     keep after symm
    //     should we check if a map is loaded (or density scoring is on) instead??
    if ( option[ OptionKeys::edensity::mapfile ].user() || user_csts_.size() > 0 ) {
      MoverOP dens( new protocols::electron_density::SetupForDensityScoringMover );
      dens->apply( pose );
    }

    // (5) initialize template history
    //     keep after symm
    TemplateHistoryOP history( new TemplateHistory(pose) );
    history->setall( initial_template_index );
    pose.data().set( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY, history );

    // (6) allowed movement
    utility::vector1<bool> allowed_to_move;
    allowed_to_move.resize(pose.size(),true);
    for ( int i=1; i<=(int)residue_sample_template_.size(); ++i ) allowed_to_move[i] = allowed_to_move[i] && residue_sample_template_[i];
    for ( int i=1; i<=(int)residue_sample_abinitio_.size(); ++i ) allowed_to_move[i] = allowed_to_move[i] && residue_sample_abinitio_[i];

    // (7) steal crystal parameters (if specified)
    if ( templates_[initial_template_index]->pdb_info() ) {
      core::io::CrystInfo ci = templates_[initial_template_index]->pdb_info()->crystinfo();
      if ( ci.A()*ci.B()*ci.C() != 0 ) {
	pose.pdb_info()->set_crystinfo( ci );
      }
    }


    // STAGE 1
    //fpd constraints are handled a little bit weird
    //  * foldtree hybridize sets chainbreak variants then applies constraints (so c-beta csts are treated properly)
    //  * after chainbreak variants are removed, csts are removed and reapplied in this function
    //  * finally after switch to fullatom CSTs are reapplied again (optionally using a different CST file)
    // If "AUTO" is specified for cen_csts automatic centroid csts are generated
    // If "AUTO" is specified for fa_csts automatic fa csts are generated
    // If "NONE" is specified for fa_csts, cen csts are used (AUTO or otherwise)
    // An empty string is treated as equivalent to "NONE"

    // fold tree hybridize
    //if ( numeric::random::rg().uniform() < stage1_probability_ ) {
    if ( true ) {

      for ( core::Size repeatstage1=0; repeatstage1 < jump_move_repeat_; ++repeatstage1 ) {

	std::string cst_fn = template_cst_fn_[initial_template_index];

	FoldTreeHybridizeOP ft_hybridize(
					 new FoldTreeHybridize(

							       initial_template_index, templates_aln_, template_weights_,
							       template_chunks_, frags_small, frags_big) ) ;

	ft_hybridize->set_constraint_file( cst_fn );
	ft_hybridize->set_constraint( cen_cst_in_ );
	ft_hybridize->set_scorefunction( stage1_scorefxn_ );

	// options
	ft_hybridize->set_increase_cycles( 1.0 );
	ft_hybridize->set_stage1_1_cycles( stage1_1_cycles_ );
	ft_hybridize->set_stage1_2_cycles( stage1_2_cycles_ );
	ft_hybridize->set_stage1_3_cycles( stage1_3_cycles_ );
        ft_hybridize->set_stage1_4_cycles( stage1_4_cycles_ );

	//ft_hybridize->set_increase_cycles( 1 );
	//ft_hybridize->set_stage1_1_cycles( 0);
	//ft_hybridize->set_stage1_2_cycles(20);
	//ft_hybridize->set_stage1_3_cycles(0);
	//ft_hybridize->set_stage1_4_cycles(0);
	ft_hybridize->set_add_non_init_chunks( add_non_init_chunks_ );
	//ft_hybridize->set_add_hetatm( add_hetatm_, hetatm_self_cst_weight_, hetatm_prot_cst_weight_ );
	ft_hybridize->set_add_hetatm(false, 0.0, 0.0);
	ft_hybridize->set_frag_1mer_insertion_weight( frag_1mer_insertion_weight_ );
	ft_hybridize->set_small_frag_insertion_weight( small_frag_insertion_weight_ );
	ft_hybridize->set_big_frag_insertion_weight( big_frag_insertion_weight_ );
	ft_hybridize->set_chunk_insertion_weight( chunk_insertion_weight_ );
	ft_hybridize->set_frag_weight_aligned( frag_weight_aligned_ );
	ft_hybridize->set_auto_frag_insertion_weight( auto_frag_insertion_weight_ );
	ft_hybridize->set_max_registry_shift( max_registry_shift_ );

	// strand pairings
	ft_hybridize->set_pairings_file( pairings_file_ );
	ft_hybridize->set_sheets( sheets_ );
	ft_hybridize->set_random_sheets( random_sheets_ );
	ft_hybridize->set_filter_templates( filter_templates_ );

	// allowed movement
	ft_hybridize->set_per_residue_controls( residue_sample_template_, residue_sample_abinitio_ );
	ft_hybridize->set_minimize_at_end( min_after_stage1_ );
	ft_hybridize->set_minimize_sf( stage2_scorefxn_ );

	// max insertion
	ft_hybridize->set_max_insertion( max_contig_insertion_ );

	// other cst stuff
	ft_hybridize->set_user_csts( user_csts_ );

	// finally run stage 1
	//std::cout << "START stage 1 with " << (*stage1_scorefxn_ )(pose) << std::endl;

	ft_hybridize->apply(pose);
	//std::cout << "END stage 1 with " << (*stage1_scorefxn_ )(pose) << std::endl;
	// get strand pairings if they exist for constraints in later stages
	strand_pairs_ = ft_hybridize->get_strand_pairs();



	//jump perturbation and minimization
	//fpd!  these docking moves should be incorporated as part of stage 1!
	if ( false ) {
	  do_intrastage_docking( pose );
	}
      } //end of repeatstage1

    } else {
      // just do frag insertion in unaligned regions
      core::pose::PoseOP chosen_templ = templates_aln_[initial_template_index];
      protocols::loops::Loops chosen_contigs = template_contigs_[initial_template_index];
      initialize_and_sample_loops(pose, chosen_templ, chosen_contigs, stage1_scorefxn_);
    }

    // realign domains to the output of stage 1
    if ( jump_move_ || realign_domains_stage2_ ) {
      TR << "Realigning template domains to stage1 pose." << std::endl;
      core::pose::PoseOP stage1pose( new core::pose::Pose( pose ) );
      align_templates_by_domain(stage1pose); // don't use stored domains; recompute parse from model
    }

    //write gdtmm to output
    if ( native_ && native_->size() ) {
      gdtmm = get_gdtmm(*native_, pose, aln_);
      core::pose::setPoseExtraScore( pose, "GDTMM_after_stage1", gdtmm);
      TR << "GDTMM_after_stage1" << F(8,3,gdtmm) << std::endl;
    }


    // 	// STAGE 2
    if (current_stage == complete) {
      // apply constraints
      if ( stage2_scorefxn_->get_weight( core::scoring::atom_pair_constraint ) != 0 ) {
	std::string cst_fn = template_cst_fn_[initial_template_index];
	if ( !keep_pose_constraint_ ) {
	  if ( ! cen_cst_in_.empty() ) {
	    setup_constraints(pose, cen_cst_in_);
	  } else {
	    setup_centroid_constraints( pose, templates_aln_, template_weights_, cst_fn );
	  }
	}
	if ( add_hetatm_ ) {
	  add_non_protein_cst(pose, *templates_aln_[initial_template_index], hetatm_self_cst_weight_, hetatm_prot_cst_weight_);
	}
	if ( strand_pairs_.size() ) {
	  add_strand_pairs_cst(pose, strand_pairs_);
	}
      }

      // add torsion constraints derived from fragments
      if ( csts_from_frags_ ) {
	if ( stage2_scorefxn_->get_weight( core::scoring::dihedral_constraint ) != 0 ) {
	  add_fragment_csts( pose );
	} else {
	  TR << "Warning! csts_from_frags is on but dihedral_constraint weight=0.  Ignoring!" << std::endl;
	}
      }

      if ( stage2_scorefxn_->get_weight( core::scoring::coordinate_constraint ) != 0 ) {
	if ( user_csts_.size() > 0 ) {
	  setup_user_coordinate_constraints(pose,user_csts_);
	}
      }

      if ( !option[cm::hybridize::skip_stage2]() ) {
	core::scoring::ScoreFunctionOP stage2_scorefxn_clone = stage2_scorefxn_->clone();

	core::scoring::methods::EnergyMethodOptions lowres_options(stage2_scorefxn_clone->energy_method_options());
	lowres_options.set_cartesian_bonded_linear(true);
	stage2_scorefxn_clone->set_energy_method_options(lowres_options);

	CartesianHybridizeOP cart_hybridize(
					    new CartesianHybridize( templates_aln_, template_weights_, template_chunks_,template_contigs_, frags_big ) );

	// default scorefunctions (cenrot-compatable)
	if ( stage2_scorefxn_!=nullptr ) cart_hybridize->set_scorefunction( stage2_scorefxn_ );
	if ( stage2pack_scorefxn_!=nullptr ) cart_hybridize->set_pack_scorefunction( stage2pack_scorefxn_ );
	if ( stage2min_scorefxn_!=nullptr ) cart_hybridize->set_min_scorefunction( stage2min_scorefxn_ );

	// set options
	cart_hybridize->set_increase_cycles( stage2_increase_cycles_ );
	cart_hybridize->set_no_global_frame( no_global_frame_ );
	cart_hybridize->set_linmin_only( linmin_only_ );
	cart_hybridize->set_seqfrags_only( seqfrags_only_ );
	cart_hybridize->set_cartfrag_overlap( cartfrag_overlap_ );
	cart_hybridize->set_skip_long_min( skip_long_min_ );
	cart_hybridize->set_cenrot( cenrot_ );
	cart_hybridize->set_fragment_probs(fragprob_stage2_, randfragprob_stage2_);

	// max insertion
	cart_hybridize->set_max_insertion( max_contig_insertion_ );

	// per-residue controls
	cart_hybridize->set_per_residue_controls( residue_sample_template_, residue_sample_abinitio_ );

	// finally run stage 2
	cart_hybridize->apply(pose);
      }

      //write gdtmm to output
      if ( native_ && native_->size() ) {
	gdtmm = get_gdtmm(*native_, pose, aln_);
	core::pose::setPoseExtraScore( pose, "GDTMM_after_stage2", gdtmm);
	TR << "GDTMM_after_stage2" << ObjexxFCL::format::F(8,3,gdtmm) << std::endl;
      }
      // get fragment history
      runtime_assert( pose.data().has( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY ) );
      history = utility::pointer::static_pointer_cast< TemplateHistory >( pose.data().get_ptr( CacheableDataType::TEMPLATE_HYBRIDIZATION_HISTORY ) );

      TR << "History :";
      for ( Size i=1; i<= history->size(); ++i ) { TR << I(4,i); }
      TR << std::endl;
      TR << "History :";
      for ( Size i=1; i<= history->size(); ++i ) { TR << I(4, history->get(i)); }
      TR << std::endl;

      core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
      mm->set_bb  ( true );
      mm->set_chi ( true );
      mm->set_jump( true );

      //fpd -- in positions where no fragment insertions are allowed, also allow no minimization
      if ( residue_sample_template_.size() > 0 && residue_sample_abinitio_.size() > 0 ) {
	for ( int i=1; i<=(int)nres_tgt; ++i ) {
	  if ( !residue_sample_template_[i] && !residue_sample_abinitio_[i] ) {
	    TR.Trace << "locking residue " << i << std::endl;
	    mm->set_bb  ( i, false );
	    mm->set_chi ( i, false );
	  }
	}
      }

      if ( core::pose::symmetry::is_symmetric(pose) ) {
	core::pose::symmetry::make_symmetric_movemap( pose, *mm );
      }

      // stage "2.5" .. minimize with centroid energy + full-strength cart bonded
      if ( !option[cm::hybridize::skip_stage2]() ) {
	core::optimization::MinimizerOptions options_lbfgs( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
	core::optimization::CartesianMinimizer minimizer;
	auto n_min_cycles =(Size) (200.*stage25_increase_cycles_);
	options_lbfgs.max_iter(n_min_cycles);
	(*stage2_scorefxn_)(pose); minimizer.run( pose, *mm, *stage2_scorefxn_, options_lbfgs );
      }

      // STAGE 3: RELAX
      if ( batch_relax_ > 0 ) {
	// set disulfides before going to FA
	if ( disulf_file_.length() > 0 ) {
	  // delete current disulfides
	  for ( int i=1; i<=(int)pose.size(); ++i ) {
	    if ( !pose.residue(i).is_protein() ) continue;
	    if ( pose.residue(i).type().is_disulfide_bonded() ) {
	      core::conformation::change_cys_state( i, "CYS", pose.conformation() );
	    }
	  }

	  // manually add new ones
	  TR << " add disulfide: " << std::endl;
	  basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].value(disulf_file_);
	  core::pose::initialize_disulfide_bonds(pose);
	  // must reset this since initialize_disulfide_bonds is used in pose io
	  basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].deactivate();
	  basic::options::option[ basic::options::OptionKeys::in::fix_disulf ].to_default(); // reset to the default value
	} else {
	  pose.conformation().detect_disulfides();
	}

	protocols::moves::MoverOP tofa( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::FA_STANDARD ) );
	tofa->apply(pose);

	// apply fa constraints
	std::string cst_fn = template_cst_fn_[initial_template_index];
	if ( fa_scorefxn_->get_weight( core::scoring::atom_pair_constraint ) != 0 ) {
	  if ( !keep_pose_constraint_ ) {
	    if ( ! fa_cst_in_.empty() ) {
	      setup_constraints(pose, fa_cst_in_);
	    } else {
	      setup_fullatom_constraints( pose, templates_aln_, template_weights_, cst_fn, fa_cst_fn_ );
	    }
	  } else {
	    pose.constraint_set(save_pose_constraint_set);
	  }
	  if ( add_hetatm_ ) {
	    add_non_protein_cst(pose, *templates_aln_[initial_template_index], hetatm_self_cst_weight_, hetatm_prot_cst_weight_);
	  }
	  if ( strand_pairs_.size() ) {
	    add_strand_pairs_cst(pose, strand_pairs_);
	  }
	}

	// add torsion constraints derived from fragments
	if ( csts_from_frags_ ) {
	  if ( fa_scorefxn_->get_weight( core::scoring::dihedral_constraint ) != 0 ) {
	    add_fragment_csts( pose );
	  } else {
	    TR << "Warning! csts_from_frags is on but dihedral_constraint weight=0.  Ignoring!" << std::endl;
	  }
	}

	if ( stage2_scorefxn_->get_weight( core::scoring::coordinate_constraint ) != 0 ) {
	  // note that CSTs are updated based on stage 1 movement
	  //   this may or may not be desirable
	  if ( user_csts_.size() > 0 ) {
	    setup_user_coordinate_constraints(pose,user_csts_);
	  }
	}

	if ( batch_relax_ == 1 ) {
	  // standard relax
	  // do relax _without_ ramping down coordinate constraints
	  TR << " batch_relax 1 : " << std::endl;
	  protocols::relax::FastRelax relax_prot( fa_scorefxn_, relax_repeats_ ,"NO CST RAMPING" );
	  relax_prot.min_type("lbfgs_armijo_nonmonotone");
	  relax_prot.apply(pose);
	} else {
	  // batch relax
	  // add stucture to queue
	  SilentFileOptions opts;
	  SilentStructOP new_struct = SilentStructFactory::get_instance()->get_silent_struct("binary",opts);
	  new_struct->fill_struct( pose );
	  new_struct->energies_from_pose( pose );
	  post_centroid_structs.push_back( new_struct );

	  if ( post_centroid_structs.size() == batch_relax_ ) {
	    protocols::relax::FastRelax relax_prot( fa_scorefxn_ );
	    relax_prot.min_type("lbfgs_armijo_nonmonotone");
	    relax_prot.set_force_nonideal(true);
	    relax_prot.set_script_to_batchrelax_default( relax_repeats_ );

	    // need to use a packer task factory to handle poses with different disulfide patterning
	    core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory );
	    tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
	    tf->push_back( TaskOperationCOP( new core::pack::task::operation::IncludeCurrent ) );
	    tf->push_back( TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );
	    relax_prot.set_task_factory( tf );

	    // notice! this assumes all poses in a set have the same constraints!
	    relax_prot.batch_apply(post_centroid_structs, pose.constraint_set()->clone());

	    // reinflate pose
	    post_centroid_structs[0]->fill_pose( pose );
	  } else {
	    protocols::moves::MoverOP tocen( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );
	    tocen->apply(pose);
	    need_more_samples = true;
	  }
	}
      }
      //   } else {
      // 	// no fullatom sampling
      // 	(*stage2_scorefxn_)(pose);
      //   }
      // }
      if ( native_ && native_->size() ) {
	gdtmm = get_gdtmm(*native_, pose, aln_);
	core::pose::setPoseExtraScore( pose, "GDTMM_final", gdtmm);
	TR << "GDTMM_final" << ObjexxFCL::format::F(8,3,gdtmm) << std::endl;
      }
    }
}
