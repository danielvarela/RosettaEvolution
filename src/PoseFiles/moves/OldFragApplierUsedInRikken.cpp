
#include "OldFragApplierUsedInRikken.hpp"
//#include "FragApplier.hpp"

// //keep first
// #include <core/scoring/cryst/PhenixInterface.hh>


// #include <utility/pointer/owning_ptr.fwd.hh>
// #include <utility/pointer/owning_ptr.functions.hh>

// #include <utility/pointer/boost/ReferenceCount.hh>


// #include <protocols/hybridization/HybridizeProtocolCreator.hh>
// #include <protocols/hybridization/CartesianSampler.hh>
// #include <protocols/hybridization/FragmentBiasAssigner.hh>
// #include <protocols/hybridization/util.hh>
// #include <core/pose/datacache/CacheableDataType.hh>
// #include <basic/datacache/BasicDataCache.hh>

// #include <protocols/moves/Mover.fwd.hh>
// #include <protocols/moves/Mover.hh>
// #include <protocols/moves/MoverContainer.hh>
// #include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
// #include <protocols/relax/util.hh>

// #include <protocols/simple_moves/ReturnSidechainMover.hh>
// #include <protocols/simple_moves/FragmentMover.hh>
// //#include <protocols/simple_moves/symmetry/SetupNCSMover.hh>

// #include <protocols/moves/MonteCarlo.hh>
// #include <protocols/comparative_modeling/coord_util.hh>
// #include <basic/datacache/DataMap.hh>

// #include <protocols/loops/loops_main.hh>
// #include <protocols/loops/Loops.hh>

// #include <protocols/jd2/util.hh>

// #include <core/pack/task/operation/TaskOperations.hh>
// #include <core/pack/task/TaskFactory.hh>
// // #include <protocols/simple_moves/PackRotamersMover.hh>
// // #include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
// // #include <protocols/simple_moves/RotamerTrialsMover.hh>

// //#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
// #include <core/scoring/ScoreFunction.hh>
// #include <core/scoring/constraints/ConstraintSet.hh>
// #include <core/scoring/rms_util.hh>
// #include <core/scoring/electron_density/ElectronDensity.hh>
// #include <core/scoring/Energies.hh>
// #include <core/sequence/util.hh>
// #include <core/sequence/Sequence.hh>
// #include <core/id/SequenceMapping.hh>
// #include <core/sequence/SequenceAlignment.hh>

// #include <core/pose/Pose.hh>
// #include <core/pose/PDBInfo.hh>
// #include <core/io/Remarks.hh>
// #include <core/pose/selection.hh>
// #include <core/chemical/ChemicalManager.hh>
// #include <core/import_pose/import_pose.hh>

// // dump intermediate pdb
// #include <core/io/util.hh>
// #include <core/io/pdb/pdb_writer.hh>
// #include <utility/io/ozstream.hh>

// #include <core/pose/util.hh>

// #include <core/kinematics/FoldTree.hh>
// #include <core/kinematics/MoveMap.hh>

// #include <core/scoring/ScoreFunction.hh>
// #include <core/scoring/symmetry/SymmetricScoreFunction.hh>
// #include <core/scoring/ScoreFunctionFactory.hh>
// #include <core/scoring/constraints/CoordinateConstraint.hh>
// #include <core/scoring/func/HarmonicFunc.hh>
// #include <core/scoring/constraints/AtomPairConstraint.hh>
// #include <core/scoring/func/ScalarWeightedFunc.hh>
// #include <core/scoring/func/SOGFunc.hh>

// #include <core/optimization/AtomTreeMinimizer.hh>
// #include <core/optimization/MinimizerOptions.hh>
// #include <core/kinematics/MoveMap.hh>
// #include <core/optimization/CartesianMinimizer.hh>
// #include <core/optimization/MinimizerOptions.hh>

// #include <core/scoring/ScoreFunction.hh>
// #include <core/scoring/ScoreFunctionFactory.hh>
// #include <core/scoring/methods/EnergyMethodOptions.hh>
// #include <core/scoring/Energies.hh>
// #include <core/scoring/constraints/AtomPairConstraint.hh>
// #include <core/scoring/func/USOGFunc.hh>

// #include <core/pack/task/PackerTask.hh>
// #include <core/pack/task/TaskFactory.hh>
// #include <core/pack/task/operation/TaskOperations.hh>
// #include <core/conformation/util.hh>

// #include <core/fragment/ConstantLengthFragSet.hh>
// #include <core/fragment/FragmentIO.hh>
// #include <core/fragment/FragSet.hh>
// #include <core/fragment/Frame.hh>
// #include <core/fragment/FrameIterator.hh>
// #include <core/fragment/FragData.hh>

// // symmetry
// #include <core/pose/symmetry/util.hh>
// #include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
// #include <core/conformation/symmetry/SymmetricConformation.hh>
// #include <core/conformation/symmetry/SymmetryInfo.hh>


// #include <utility/excn/Exceptions.hh>
// #include <utility/file/file_sys_util.hh>

// #include <numeric/random/random.hh>
// #include <numeric/xyzVector.hh>
// #include <numeric/xyz.functions.hh>
// #include <numeric/model_quality/rms.hh>
// #include <numeric/random/WeightedSampler.hh>

// #include <basic/options/option.hh>
// #include <basic/options/option_macros.hh>
// #include <basic/options/keys/OptionKeys.hh>
// #include <basic/options/keys/in.OptionKeys.gen.hh>
// #include <basic/options/keys/out.OptionKeys.gen.hh>
// #include <basic/options/keys/cm.OptionKeys.gen.hh>

// #include <utility/tag/Tag.hh>
// #include <utility/string_util.hh>
// #include <basic/Tracer.hh>

// #include <boost/unordered/unordered_map.hpp>
// // XSD XRW Includes
// #include <utility/tag/XMLSchemaGeneration.hh>
// #include <protocols/moves/mover_schemas.hh>


// #include <memory>

// using namespace core;
// using namespace core::scoring;
// using namespace core::scoring::constraints;


//   // transform fragment
//   void
//   FragApplier::apply_transform(
// 		  core::pose::Pose &frag,
// 		  core::Vector const &preT,
// 		  core::Vector const &postT,
// 		  numeric::xyzMatrix< core::Real > const &R
// 		  ){
//     // apply rotation to ALL atoms
//     // x_i' <- = R*x_i + com1;
//     for ( core::Size i = 1; i <= frag.size(); ++i ) {
//       for ( core::Size j = 1; j <= frag.residue_type(i).natoms(); ++j ) {
// 	core::id::AtomID id( j, i );
// 	frag.set_xyz( id, R * ( frag.xyz(id) - preT) + postT );
//       }
//     }
//   }


//   // apply endpoint csts from pose to frag
//   void
//   FragApplier::apply_fragcsts(
// 		 core::pose::Pose &working_frag,
// 		 core::pose::Pose const &pose,
// 		 core::Size startpos
// 		 ){
//     // using namespace core::scoring::constraints;

//     // working_frag.remove_constraints();
//     // int len = working_frag.size();
//     // if ( working_frag.residue(len).aa() == core::chemical::aa_vrt ) len--;

//     // bool freeze_endpoints_ = true;
//     // bool nterm = ! freeze_endpoints_ && ( (startpos == 1) || pose.fold_tree().is_cutpoint(startpos-1) );
//     // bool cterm = ! freeze_endpoints_ && ( pose.fold_tree( ).is_cutpoint(startpos+len-1) );

//     // if ( !nterm ) {
//     //   for ( core::uint j = 0; j < overlap_; ++j ) {
//     // 	for ( core::uint i = 1; i <= 3; ++i ) {
//     // 	  core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
//     // 	  working_frag.add_constraint(
//     // 				      scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new CoordinateConstraint(
//     // 																	core::id::AtomID(i,j+1),
//     // 																	core::id::AtomID(2,working_frag.size()),
//     // 																	pose.residue(startpos+j).atom(i).xyz(),
//     // 																	fx
//     // 																	) ) )
//     // 				      );
//     // 	}
//     //   }
//     // }
//     // if ( !cterm ) {
//     //   for ( int j=len-overlap_; j<len; ++j ) {
//     // 	for ( int i=1; i<=3; ++i ) {
//     // 	  core::scoring::func::FuncOP fx( new core::scoring::func::HarmonicFunc( 0.0, 1.0 ) );
//     // 	  working_frag.add_constraint(
//     // 				      scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new CoordinateConstraint(
//     // 																	core::id::AtomID(i,j+1),
//     // 																	core::id::AtomID(2,working_frag.size()),
//     // 																	pose.residue(startpos+j).atom(i).xyz(),
//     // 																	fx
//     // 																	) ) )
//     // 				      );
//     // 	}
//     //   }
//     // }
//   }


// // get frag->pose transform, return RMS
// // default overlap_ = 2
// core::Real
// FragApplier::get_transform(
// 	core::pose::Pose const &pose,
// 	core::pose::Pose const &frag,
// 	core::Size startpos,
// 	core::Vector &preT,
// 	core::Vector &postT,
// 	numeric::xyzMatrix< core::Real > &R
// ){
// 	int aln_len = overlap_;
// 	int len = frag.size();
// 	if ( frag.residue(len).aa() == core::chemical::aa_vrt ) len--;
// 	ObjexxFCL::FArray1D< core::Real > ww( 2*4*aln_len, 1.0 );
// 	ObjexxFCL::FArray2D< core::Real > uu( 3, 3, 0.0 ), init_coords( 3, 2*4*aln_len ), final_coords( 3, 2*4*aln_len );
// 	preT = numeric::xyzVector< core::Real >(0,0,0);
// 	postT = numeric::xyzVector< core::Real >(0,0,0);

//     bool freeze_endpoints_ = true;
// 	bool nterm = ! freeze_endpoints_ && ( (startpos == 1) || pose.fold_tree().is_cutpoint(startpos-1) );
// 	bool cterm = ! freeze_endpoints_ && ( pose.fold_tree( ).is_cutpoint(startpos+len-1) );

// 	// grab coords from input pose
// 	for ( int ii=-aln_len; ii<aln_len; ++ii ) {
// 		int i = (ii>=0) ? (nterm?len-ii-1:ii) : (cterm?-ii-1:len+ii);
// 		numeric::xyzVector< core::Real > x_1 = pose.residue(startpos+i).atom(1).xyz();
// 		numeric::xyzVector< core::Real > x_2 = pose.residue(startpos+i).atom(2).xyz();
// 		numeric::xyzVector< core::Real > x_3 = pose.residue(startpos+i).atom(3).xyz();
// 		numeric::xyzVector< core::Real > x_4 = pose.residue(startpos+i).atom(4).xyz();
// 		postT += x_1+x_2+x_3+x_4;
// 		for ( int j=0; j<3; ++j ) {
// 			init_coords(j+1,4*(ii+aln_len)+1) = x_1[j];
// 			init_coords(j+1,4*(ii+aln_len)+2) = x_2[j];
// 			init_coords(j+1,4*(ii+aln_len)+3) = x_3[j];
// 			init_coords(j+1,4*(ii+aln_len)+4) = x_4[j];
// 		}
// 	}
// 	postT /= 2.0*4.0*aln_len;
// 	for ( int ii=0; ii<2*4*aln_len; ++ii ) {
// 		for ( int j=0; j<3; ++j ) init_coords(j+1,ii+1) -= postT[j];
// 	}

// 	// grab new coords from frag
// 	for ( int ii=-aln_len; ii<aln_len; ++ii ) {
// 		int i = (ii>=0) ? (nterm?len-ii-1:ii) : (cterm?-ii-1:len+ii);
// 		numeric::xyzVector< core::Real > x_1 = frag.residue(i+1).atom(1).xyz();
// 		numeric::xyzVector< core::Real > x_2 = frag.residue(i+1).atom(2).xyz();
// 		numeric::xyzVector< core::Real > x_3 = frag.residue(i+1).atom(3).xyz();
// 		numeric::xyzVector< core::Real > x_4 = frag.residue(i+1).atom(4).xyz();
// 		preT += x_1+x_2+x_3+x_4;
// 		for ( int j=0; j<3; ++j ) {
// 			final_coords(j+1,4*(ii+aln_len)+1) = x_1[j];
// 			final_coords(j+1,4*(ii+aln_len)+2) = x_2[j];
// 			final_coords(j+1,4*(ii+aln_len)+3) = x_3[j];
// 			final_coords(j+1,4*(ii+aln_len)+4) = x_4[j];
// 		}
// 	}
// 	preT /= 2.0*4.0*aln_len;
// 	for ( int ii=0; ii<2*4*aln_len; ++ii ) {
// 		for ( int j=0; j<3; ++j ) final_coords(j+1,ii+1) -= preT[j];
// 	}

// 	numeric::Real ctx;
// 	float rms;
// 	numeric::model_quality::findUU( final_coords, init_coords, ww, 2*4*aln_len, uu, ctx );
// 	numeric::model_quality::calc_rms_fast( rms, final_coords, init_coords, ww, 2*4*aln_len, ctx );

// 	R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
// 	R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
// 	R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );

// 	return ((core::Real)rms);
// }

//   bool
//   FragApplier::apply_frame (
// 		    core::pose::Pose & pose /*, core::fragment::Frame frame*/
// ){
// 	using core::pack::task::operation::TaskOperationCOP;

// 	//core::Size start = frame.start(),len = frame.length();
// 	core::Size start = 2,len = 9;
// 	runtime_assert( overlap_>=1 && overlap_<=len/2); // at least two-residue overlap, 4 mer at least

// 	// set the frame's insert point to 1
// 	//	frame.shift_to(1);

// 	// // see if the pose has NCS
// 	// simple_moves::symmetry::NCSResMappingOP ncs;
// 	// if ( pose.data().has( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ) {
// 	// 	ncs = ( utility::pointer::static_pointer_cast< simple_moves::symmetry::NCSResMapping > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ));
// 	// }

// 	protocols::simple_moves::SwitchResidueTypeSetMover to_fa("fa_standard");
// 	protocols::simple_moves::SwitchResidueTypeSetMover to_cen("centroid");

// 	// make subpose at this position
// 	// we assume the frame does not cross a jump
//         //Maybe this can be done during differential evolution to avoid this step??
// 	//During DE always apply to the same frame, but different FragData! (= each individual)
// 	core::pose::Pose frag = core::pose::Pose(pose, 2, 10);
// 	// for ( int i=0; i<(int)len; ++i ) frag.append_residue_by_bond( pose.residue( start+i ) );
// 	// for ( int i=0; i<(int)len; ++i ) {
// 	// 	if ( frag.residue(i+1).aa() == chemical::aa_cys && frag.residue(i+1).has_variant_type( chemical::DISULFIDE ) ) {
// 	// 		TR << "temp remove dslf at " << i+1 << std::endl;
// 	// 		conformation::change_cys_state( i+1, "CYS", frag.conformation() );
// 	// 	}
// 	// }

// 	for ( int i=0; i<(int)len; ++i ) core::conformation::idealize_position(i+1, frag.conformation());

// 	core::Vector preT, postT;
// 	numeric::xyzMatrix< core::Real > R;
// 	core::Size maxtries,frag_toget=0;
// 	core::Real rms_cutoff_ = 1.5;
// 	core::Real rms=rms_cutoff_;
//         bool bbmove_ = true;

// 	// set up minimizer
// 	core::optimization::AtomTreeMinimizer rbminimizer;
// 	core::optimization::MinimizerOptions options_rb( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
// 	options_rb.max_iter(20);
// 	core::kinematics::MoveMap mm_rb;
// 	mm_rb.set_bb ( bbmove_ );
// 	mm_rb.set_chi ( true );
// 	mm_rb.set_jump ( true );


// 	core::scoring::ScoreFunctionOP fa_scorefxn_ = core::scoring::get_score_function();
//  core::scoring::ScoreFunctionOP scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function(std::string("ref2015_cst").c_str());


// 	// set up packer
// 	core::scoring::ScoreFunctionOP nonsymm_fa_scorefxn = core::scoring::symmetry::asymmetrize_scorefunction(*fa_scorefxn_);
// 	core::scoring::ScoreFunctionOP nonsymm_cen_scorefxn = core::scoring::symmetry::asymmetrize_scorefunction(*scorefxn_);

// 	nonsymm_fa_scorefxn->set_weight( core::scoring::coordinate_constraint , 1 );

// 	core::pack::task::TaskFactoryOP main_task_factory( new core::pack::task::TaskFactory );
// 	//	boost::shared_ptr<core::pack::task::operation::TaskOperation const > restrict( );
// 	main_task_factory->push_back( boost::shared_ptr<core::pack::task::operation::TaskOperation>(new core::pack::task::operation::RestrictToRepacking)  );

// 	std::string selection_bias_ = "density";
//         bool fullatom_ = true;

//         if ( selection_bias_ == "density" ) {
// 		// scorefunctions
// 	  nonsymm_fa_scorefxn->set_weight( core::scoring::elec_dens_fast , 20 );

// 		core::scoring::ScoreFunctionOP densonly( new core::scoring::ScoreFunction() );
// 		densonly->set_weight( core::scoring::elec_dens_fast, 5.0 );

// 		// prepare fragment
// 		core::pose::addVirtualResAsRoot(frag);
// 		if ( !fullatom_ ) to_fa.apply( frag );

// 		core::Size nattempts = 0;
// 		core::Real best_dens_score = 1e30, best_rms = 1e30, best_rms_aftermin = 1e30;
// 		//for ( core::uint i = 1; i <= frame.nr_frags(); ++i ) {
// 		core::uint i = 1;
// 		core::pose::Pose working_frag = frag;

// 		//	frame.apply( i, working_frag );
// 			rms = get_transform( pose,  working_frag,  start, preT, postT, R); // rms straight from comparing frag to pose at region with nres overlap_
// 			if ( rms <= rms_cutoff_ ) {
// 				nattempts++;

// 				// orient
// 				apply_transform( working_frag, preT, postT, R);

// 				// cenmin+repack+min
// 				to_cen.apply( working_frag );
// 				(*densonly)(working_frag);
// 				rbminimizer.run( working_frag, mm_rb, *densonly, options_rb );
// 				to_fa.apply( working_frag );

// 				core::Real dens_score;
// 			        // apply csts from pose to frag at endpoints
// 				apply_fragcsts( working_frag, pose, start );
// 				(*nonsymm_fa_scorefxn)(working_frag);
// 				rbminimizer.run( working_frag, mm_rb, *nonsymm_fa_scorefxn, options_rb );
// 				dens_score = (*nonsymm_fa_scorefxn) (working_frag);

// 				if ( dens_score<best_dens_score ) {
// 					best_rms_aftermin = get_transform( pose,  working_frag,  start, preT, postT, R); // rms after minimization
// 					best_rms = rms; // this rms is before minimization, should return after minimization
// 					best_dens_score = dens_score;
// 					frag = working_frag;
// 				}
// 				//if ( nattempts >=25 ) break;  //fpd within helices often all fragments will match

// 			} // if rms <= rms_cutoff_
// 			//} // try all the frags in the position until some have rms < cutoff

// 		if ( nattempts > 0 ) {
// 		  std::cout << "Best frag ( out of "<< nattempts << ") with score="  << best_dens_score << "  rms=" << best_rms << "  rms_aftermin=" << best_rms_aftermin << std::endl;

// 		  for ( core::Size i = 0; i < len; ++i ) {
// 				bool dslf_i = ( pose.residue(start+i).aa() == chemical::aa_cys && pose.residue(start+i).has_variant_type( chemical::DISULFIDE ) );

// 				// remember dslf connectivity
// 				if ( dslf_i ) {
// 					for ( int j=1; j<=(int)pose.residue(start+i).natoms(); ++j ) {  // copy everything but HG
// 					  std::cout << "restore dslf at " << start+i << std::endl;
// 						pose.set_xyz( core::id::AtomID( j,start+i ), frag.residue(i+1).xyz( j ) );
// 					}
// 				} else {
// 					pose.replace_residue( start+i, frag.residue(i+1), false );
// 				}
// 			}
// 		} else { // best scoring fragment saved
// 		  std::cout << "No acceptable fragments" << std::endl;
// 		  //	frame.shift_to(start);
// 			return false;
// 		}
// 	}

// 	// // apply to NCS-symmetric copies
// 	// if ( ncs ) {
// 	// 	for ( core::uint j = 1; j <= ncs->ngroups(); ++j ) {
// 	// 		bool all_are_mapped = true;
// 	// 		for ( Size k= 0; k< len && all_are_mapped; ++k ) {
// 	// 			all_are_mapped &= (ncs->get_equiv( j,start+k )!=0);
// 	// 		}
// 	// 		if ( !all_are_mapped ) continue;

// 	// 		core::Size remap_start = ncs->get_equiv( j, start );
// 	// 		get_transform( pose,  frag,  remap_start, preT, postT, R);
// 	// 		apply_transform( frag, preT, postT, R);
// 	// 		for ( Size i = 0; i < len; ++i ) {
// 	// 			bool dslf_i = ( pose.residue(remap_start+i).aa() == chemical::aa_cys && pose.residue(remap_start+i).has_variant_type( chemical::DISULFIDE ) );
// 	// 			if ( dslf_i ) {
// 	// 				for ( int k=1; k<=(int)pose.residue(remap_start+i).natoms(); ++k ) {  // copy everything but HG
// 	// 					TR << "(ncs) restore dslf at " << remap_start+i << std::endl;
// 	// 					pose.set_xyz( core::id::AtomID( k,remap_start+i ), frag.residue(i+1).xyz( k ) );
// 	// 				}
// 	// 			} else {
// 	// 				pose.replace_residue( remap_start+i, frag.residue(i+1), false );
// 	// 			}
// 	// 		}
// 	// 	}
// 	// }

// 	// // restore frame
// 	// frame.shift_to(start);

//     return true;
// }
