
#include "PoseFunction.hpp"

#include <ctime>
#include <iostream>

#include "../Algorithm/DE_types.hpp"
#include "FitFunction.hpp"
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/CartesianMinimizerMap.hh>
#include <core/optimization/CartesianMultifunc.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/chemical/ResidueType.hh>


#include <core/pose/copydofs/util.hh>
#include <core/pose/copydofs/CopyDofs.hh>
#include <core/pose/copydofs/CopyDofsInfo.hh>
#include <protocols/simple_moves/CopyDofMover.hh>
#include <core/pose/copydofs/util.hh>


#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/util.hh>

#include <core/scoring/rms_util.hh>
#include <protocols/loops/Loops.hh>
#include <core/scoring/rms_util.hh>
#include <core/pose/Pose.hh>
#include <core/sequence/util.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/pose/util.hh>

#include <protocols/hybridization/DDomainParse.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <core/conformation/util.hh>
#include <core/conformation/Residue.hh>

#include <core/fragment/Frame.hh>
#include <core/fragment/BBTorsionAndAnglesSRFD.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>


double PoseImprovedFunction::score(Individual& ind) {
  // using namespace core;
  // using namespace core::fragment;
  // using namespace core::conformation;
  // using namespace fragment;
  // using core::fragment::BBTorsionSRFD;

  // Individual ind_converted = convert(ind);
  // fill_pose(pose_, ind_converted, ss);
  // // core::fragment::FragData* dat = new core::fragment::FragData( core::fragment::SingleResidueFragDataOP( new BBTorsionSRFD , 3 ) );


  // // core::fragment::SingleResidueFragDataOP sr = core::fragment::SingleResidueFragDataOP();
  // BBTorsionAndAnglesSRFD bbtorsion ;
  // // utility::pointer::shared_ptr< core::fragment::SingleResidueFragData > sr =   utility::pointer::shared_ptr< core::fragment::SingleResidueFragData>( new BBTorsionAndAnglesSRFD()  );

  // core::Size s_ = 5;


  // FrameOP frame( new Frame(  (range.j - range.i)    ) );
  // core::fragment::FragDataCOP dat = FragDataCOP(core::fragment::FragDataOP(new core::fragment::FragData(
  // 													bbtorsion.clone(), pose_->size()													   )));

  
  //   FragDataOP current_fragment( nullptr );
  //   //    current_fragment = FragDataOP( new AnnotatedFragData( ) );
  //   //SingleResidueFragDataOP sr = utility::pointer::dynamic_pointer_cast<SingleResidueFragDataOP >(bbtorsion);
  //   //   core::fragment::FragData frag( bbtorsion.clone() );

  //   //current_fragment = frag.clone();
  //   current_fragment = FragDataOP( new AnnotatedFragData( std::string("123j"), range.i ) );

  //   for (int i = range.i; i < range.j; ++i) {
  //     //      Residue r   =  pose_->get_residue(i);
  //     std::string pdbid = "123j";
  //     core::Size index        = i;
  //     char aa           = pose_->aa(i);
  //     char ss           = pose_->secstruct()[i - 1];
  //     core::Real phi          = pose_->phi(i);
  //     core::Real psi          = pose_->psi(i);
  //     core::Real omega        = pose_->omega(i);

  //     utility::pointer::shared_ptr< BBTorsionSRFD > res_torsions( new BBTorsionSRFD(3,ss,aa) ); // 3 protein torsions
  //     res_torsions->set_torsion   ( 1, phi   ); // ugly numbers 1-3, but pose.set_phi also uses explicit numbers
  //     res_torsions->set_torsion   ( 2, psi   );
  //     res_torsions->set_torsion   ( 3, omega );
  //     res_torsions->set_secstruct ( ss );


  //     res_torsions->show(std::cout);
  //     // Add residue to fragment
  //     current_fragment->add_residue( res_torsions );
  //   if (current_fragment) {
  //      current_fragment->show(std::cout);
  //   }

  //   }

  //   int parada_;
  //   std::cout << "fin fragment" << std::endl;
  //   std::cin >> parada_;

  // // core::fragment::FragDataOP fData();
  //   //  core::fragment::FrameOP frame = core::fragment::FrameOP(new core::fragment::Frame( range.i, (range.j - range.i) ) );
  // // FrameOP frame( new Frame( range.i, fData, (range.j - range.i) ) ) ;

  // frame->steal( *pose_ );
  // frame->show_classic(std::cout);
  // frame->apply( *mm_rb, range.i,   *pose_ );

  // // //newwwww
  // // core::pose::Pose frag = core::pose::Pose(*pose_, range.i, range.j + 1);
  // // // //idealize
  // // // for ( int i=0; i<(int)D_ - 1; ++i ) core::conformation::idealize_position(i+1, frag.conformation());

  // // to_cen->apply( frag );
  // // (*densonly)(frag);
  // // rbminimizer->run( frag, *mm_rb, *densonly, *options_rb );
  // // to_fa->apply( frag );
  // // pack_mover->apply( frag );

  // // int r = 1;
  // // for ( core::Size i = range.i; i < range.j; ++i ) {
  // //   pose_->replace_residue( i, frag.residue(r) , true );
  // //   r++;
  // // }

  // double result = (*scorefxn)(*pose_);
  // ind[D()] = result;
  // return result;
  return 1;
}



double
PoseFragmentFunction::run_frag_mover(core::pose::Pose& pose) {
  frag_mover->apply(pose);
  double result = SCORE_ERROR_FIXED +  (*scorefxn)(pose);
  return result;
}

double
PoseFragmentFunction::apply_fragment_stage(Individual& ind) {
  stats.increment("de_individuals");
  core::pose::PoseOP inner_pose = native->clone();
  fill_pose(inner_pose, ind, ss);
  double result = SCORE_ERROR_FIXED + (*scorefxn)(*inner_pose);
  double compare_antes = (-1 * SCORE_ERROR_FIXED) + result;
  //    std::cout << "antes : " <<  (-1 * SCORE_ERROR_FIXED) + result  << " ";
  stats.record("avg_before_de", (-1 * SCORE_ERROR_FIXED) + result );
  ind.score = result;
  core::pose::PoseOP pose_backup = inner_pose->clone();
  double result_after_frags = run_frag_mover(*inner_pose);
      
  //std::cout << "despues : " <<  (-1 * SCORE_ERROR_FIXED) + result_after_frags  << std::endl;
  double compare_result =  (-1 * SCORE_ERROR_FIXED) + result_after_frags ;
  //std::cout << "compare result " << compare_result << " " << compare_antes << " -> " << (compare_result < result) << std::endl;
  //if ( std::abs(result_after_frags - result) > 0.00001) {
  if ( compare_result < compare_antes) {
    PoseScoreFunction::pose_to_ind(inner_pose, ind);
    ind.score = result_after_frags;
    stats.increment("improved_after_de");
    stats.record("improved_amt_de",  std::abs( result_after_frags - result ) );
    result = result_after_frags;
    //std::cout << "final: " <<  (-1 * SCORE_ERROR_FIXED) + ind.score  << std::endl;
  } else {
    pose_to_ind(pose_backup, ind);
    ind.score = result;
  }

  stats.record("avg_after_de", (-1 * SCORE_ERROR_FIXED) +  ind.score );

  stats.record("total_tries_de", stats.obtain("total_tries") );
  stats.erase("total_tries");

  stats.record("accepted_fragments_de", stats.obtain("accepted_fragments") );
  stats.erase("accepted_fragments");

  return result;
}

double
PoseDensityFunction::score(Individual& ind) {
  core::pose::PoseOP inner_pose = native->clone();
  PoseDensityFunction::fill_pose(inner_pose, ind, ss);
  // to_all->apply(*inner_pose);
  // protocols::relax::FastRelax my_relax( scorefxn, 1 );
  // core::kinematics::MoveMapOP mm;
  // mm = core::kinematics::MoveMapOP( new core::kinematics::MoveMap);
  // mm->set_bb ( false );
  // mm->set_chi ( true );
  // mm->set_jump ( false );
  // my_relax.set_movemap(mm);
  // my_relax.apply(*inner_pose);
  dockindens->apply(*inner_pose);
  double result = SCORE_ERROR_FIXED + (*scorefxn)(*inner_pose);
  ind.score = result;
  return result;
}

double
PoseDensityFunction::run_frag_mover(core::pose::Pose& pose) {
  frag_mover->apply(pose);

  //to_all->apply(pose);
  // protocols::relax::FastRelax my_relax( scorefxn , 1);
  // core::kinematics::MoveMapOP mm;
  // mm = core::kinematics::MoveMapOP( new core::kinematics::MoveMap);
  // mm->set_bb ( false );
  // mm->set_chi ( true );
  // mm->set_jump ( false );
  // my_relax.set_movemap(mm);
  // my_relax.apply(pose);

  dockindens->apply(pose );
  return SCORE_ERROR_FIXED +  (*scorefxn)(pose);
}


void
PoseDensityFunction::fill_pose(core::pose::PoseOP& p, const Individual& ind, std::string ss_) {
  // from vars to pose
  core::pose::copydofs::CopyDofsInfo dofs_info_container = copy_dofs->copy_dofs_info( *native);
  dofs_info_container.clear();
  if (vector_info_native.size() != ind.vars.size()) {
    std::cout << "size diff " << vector_info_native.size() << " != " << ind.vars.size() << std::endl;
  }
  for (int i = 0; i < ind.vars.size(); i++) {
    dofs_info_container.push_back(std::pair< core::id::DOF_ID, core::Real >( vector_info_native[i+1].first, ind.vars[i]));
  }

  copy_dofs->set_copy_dofs_info(*native, dofs_info_container);
  copy_dofs->apply(*p);
}

void
PoseDensityFragmentFunction::fill_pose(core::pose::PoseOP& p, const Individual& ind, std::string ss_) {
  // from vars to pose
  core::pose::copydofs::CopyDofsInfo dofs_info_container = copy_dofs->copy_dofs_info( *native);
  dofs_info_container.clear();
  if (vector_info_native.size() != ind.vars.size()) {
    std::cout << "size diff " << vector_info_native.size() << " != " << ind.vars.size() << std::endl;
  }
  for (int i = 0; i < ind.vars.size(); i++) {
    dofs_info_container.push_back(std::pair< core::id::DOF_ID, core::Real >( vector_info_native[i+1].first, ind.vars[i]));
  }

  core::pose::copydofs::apply_dofs(*p, dofs_info_container);
  // copy_dofs->set_copy_dofs_info(*p, dofs_info_container);
  // copy_dofs->apply(*p);
  //  PoseFragmentFunction::fill_pose(p, ind, ss_);
}

void
PoseDensityFragmentFunction::apply_relax_stage(core::pose::PoseOP & inner_pose) {
  to_all->apply(*inner_pose);
  protocols::relax::FastRelax my_relax( scorefxn_fa , 1);
  my_relax.cartesian(true);
  core::kinematics::MoveMapOP mm;
  mm = core::kinematics::MoveMapOP( new core::kinematics::MoveMap);
  mm->set_bb ( true );
  mm->set_chi ( true );
  mm->set_jump ( false );
  my_relax.set_movemap(mm);
  my_relax.apply(*inner_pose);
  to_cen->apply(*inner_pose);
  dockindens->apply(*inner_pose);
}

double
PoseDensityFragmentFunction::apply_fragment_stage(Individual& ind) {
  //scorefxn->set_weight( core::scoring::elec_dens_fast, 30.0 );
  stats.increment("de_individuals");
  double begin_score = (-1 * SCORE_ERROR_FIXED) +  ind.score;
  core::pose::PoseOP pose_with_frags = native->clone();
  core::pose::PoseOP inner_pose = native->clone();
  //fill_pose(inner_pose, ind, ss);
  core::pose::copydofs::CopyDofsInfo dofs_info_container_in = copy_dofs->copy_dofs_info( *native);
  dofs_info_container_in.clear();
  // if (vector_info_native.size() != ind.vars.size()) {
  //   std::cout << "size diff " << vector_info_native.size() << " != " << ind.vars.size() << std::endl;
  // }
  for (int i = 0; i < ind.vars.size(); i++) {
    dofs_info_container_in.push_back(std::pair< core::id::DOF_ID, core::Real >( vector_info_native[i+1].first, ind.vars[i]));
  }
  core::pose::copydofs::apply_dofs(*inner_pose, dofs_info_container_in);

 //copy_dofs->set_copy_dofs_info(*inner_pose, dofs_info_container_in);
  //copy_dofs->apply(*inner_pose);

  int size_pose = inner_pose->total_residue();
  //boost::shared_ptr<HybridRosettaInserter> HybridMover = boost::dynamic_pointer_cast<HybridRosettaInserter>(frag_mover);
  core::pose::PoseOP pose_backup = inner_pose->clone();

  //core::pose::Pose my_pose = core::pose::Pose(*inner_pose);
  dockindens->apply(*inner_pose);
  double result = SCORE_ERROR_FIXED + (*scorefxn)(*inner_pose);

  //  double result_1 = (*scorefxn)(my_pose);
  //  core::scoring::ScoreFunctionOP scorefxn_nodensity = core::scoring::ScoreFunctionFactory::create_score_function("score3");
  // scorefxn_nodensity->set_weight( core::scoring::elec_dens_fast, 0.0 );
  //  double result_nodensity = (*scorefxn_nodensity)(my_pose);
  //  std::cout << "result " << result_1 << " w/o density " << result_nodensity << std::endl;
  // if ( size_pose < inner_pose->total_residue()) {
  //     core::pose::PoseOP recortado = core::pose::PoseOP(new core::pose::Pose(my_pose,1, size_pose));
  //     inner_pose = recortado->clone();
  //   }
  //   if (inner_pose->is_centroid()) {
  //   } else {
  //     protocols::moves::MoverOP tocen( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );
  //     tocen->apply(my_pose);
  //   }
  // inner_pose->energies().show_totals(std::cout);
  //std::cout << std::endl;
  // core::scoring::EnergyMap energy_map = inner_pose->energies().total_energies();
  //  std::cout << "density " << std::to_string(energy_map[core::scoring::ScoreType::elec_dens_fast]) << " " << std::endl;
  double compare_antes = (-1 * SCORE_ERROR_FIXED) + result;
  //  std::cout << "FRAGS antes : " <<  (-1 * SCORE_ERROR_FIXED) + result  << " ";
  stats.record("avg_before_de", (-1 * SCORE_ERROR_FIXED) + result );
  ind.score = result;

  //return result;
  //  std::cout << "begin : " << begin_score << " score recalc : " << (-1 * SCORE_ERROR_FIXED) + result  << std::endl;;
  //  double result_after_frags = run_frag_mover(my_pose);
  dockindens->apply(*inner_pose);
  dockindens->apply(*pose_with_frags);
  // std::cout << "iner pose before " << compare_antes << std::endl;

  //  apply_relax_stage(inner_pose);
  frag_mover->apply(*inner_pose);
  //apply_relax_stage(inner_pose);
  //std::cout << "iter 1 finish" << std::endl;
  // frag_mover->apply(*inner_pose);
  // apply_relax_stage(inner_pose);
  // std::cout << "iter 2 finish" << std::endl;
  // frag_mover->apply(*inner_pose);
  // apply_relax_stage(inner_pose);
  // std::cout << "iter 3 finish" << std::endl;
  // frag_mover->apply(*inner_pose);
  // apply_relax_stage(inner_pose); 
  // std::cout << "iter 4 finish" << std::endl;
 
  //inner_pose->dump_pdb("inner_pose_after.pdb");
  //dockindens->apply(*inner_pose);
  //  exit(1);
  // inner_pose->set_phi(1, 0);
  // inner_pose->set_psi(inner_pose->total_residue() - 1, 0);
  // //  inner_pose->constraint_set(HybridMover->hybrid_mover->save_pose_constraint_set);
  //vec.push_back(pose_with_frags->total_residue());
  // MinimizerMap min_map;
  // core::kinematics::MoveMapOP mm;
  // mm = core::kinematics::MoveMapOP( new core::kinematics::MoveMap);
  // mm->set_bb ( false );
  // mm->set_chi ( true );
  // mm->set_jump ( false );
  // min_map.setup( *inner_pose, *mm );
  // scorefxn->setup_for_minimizing( *inner_pose, min_map );
  // core::optimization::CartesianMinimizerMap min_map;
  // min_map.setup( *inner_pose, *mm );
  // core::optimization::Multivec  vars( min_map.ndofs() );
  // min_map.copy_dofs_from_pose( *inner_pose, vars );
  // std::cout << "dofs size " <<  std::to_string(vars.size()) << " " << inner_pose->total_residue() << std::endl ;
  // min_map.copy_dofs_to_pose( *pose_with_frags, vars );
  // core::optimization::Multivec  vars_2( min_map.ndofs() );
  // min_map.copy_dofs_from_pose( *pose_with_frags, vars_2 );
  // core::pose::ResMap res_map_2;
  // std::map < core::id::AtomID , core::id::AtomID > atom_id_map;
  // core::pose::copydofs::setup_atom_id_map(atom_id_map, res_map, *inner_pose);
  // antes copiaba aqui vector_info, ahora despues del if compare result
  //  std::cout << "dofs size " << vector_info.size() << std::endl;;
  // for ( auto const & elem : vector_info ) {
  //   core::Real dof_value = elem.second; // Index in the little "chunk" or "scratch" pose
  //   std::cout << dof_value << " , ";
  // }
  // std::cout << std::endl;
  // core::optimization::Multivec  dEros_dvars;
  // (*scorefxn)(*inner_pose);  // score pose first
  // scorefxn->setup_for_minimizing( *inner_pose, min_map );
  // core::optimization::CartesianMultifunc f_ros( *inner_pose, min_map, *scorefxn, false, false );
  // f_ros.dfunc( vars, dEros_dvars );

  // asi funciona, pero utilizando sus funciones internas no, muy raro!
  core::pose::copydofs::copy_dofs(*pose_with_frags, *inner_pose);

  double result_after_frags = SCORE_ERROR_FIXED + (*scorefxn)(*inner_pose);
  //double result_after_frags = SCORE_ERROR_FIXED + (*scorefxn)(*pose_with_frags);
  // std::cout << "iner pose after " << (-1 * SCORE_ERROR_FIXED) + result_after_frags << std::endl;

  //double result_check = SCORE_ERROR_FIXED + (*scorefxn)(*pose_with_frags);

  // std::cout << "result check " <<  (-1 * SCORE_ERROR_FIXED) + result_check << std::endl;
  //  dockindens->apply(*pose_with_frags);
  // double result_check_dock = SCORE_ERROR_FIXED + (*scorefxn)(*pose_with_frags);
  // std::cout << "result check dock " <<  (-1 * SCORE_ERROR_FIXED) + result_check_dock << std::endl;
  // pose_with_frags->dump_pdb("check_pose_after.pdb");
  // exit(1);

  double compare_result =  (-1 * SCORE_ERROR_FIXED) + result_after_frags ;
  // std::cout << "FRAGS despues : " <<  (-1 * SCORE_ERROR_FIXED) + result_after_frags  << std::endl;
  if ( compare_result < compare_antes) {
    //if ( std::abs(result_after_frags - result) > 0.00001) {
    // std::cout << "ENTRA AQUI " << std::endl;
    // double result_2 = SCORE_ERROR_FIXED + (*scorefxn)(*pose_with_frags);
    // std::cout << "RECALCULATE (pose) : " <<  (-1 * SCORE_ERROR_FIXED) + result_2  << std::endl;
    //PoseScoreFunction::pose_to_ind(pose_with_frags, ind);
    std::map<core::Size, core::Size > res_map;
    for (int i = 1; i <= pose_with_frags->total_residue() - 1; i++) {
      res_map[i] = i;
    }
    protocols::simple_moves::CopyDofMover copy_dofs_in( *inner_pose, res_map);
    copy_dofs_in.apply( *inner_pose );
    core::pose::copydofs::CopyDofsInfo dofs_info_container = copy_dofs_in.copy_dofs_info( *inner_pose );
    utility::vector1< std::pair< core::id::DOF_ID, core::Real > > vector_info =  dofs_info_container.dofs_info();

    std::vector<double> dofs_vars;
    for ( auto const & elem : vector_info ) {
      core::Real dof_value = elem.second; // Index in the little "chunk" or "scratch" pose
      dofs_vars.push_back(dof_value);
    }
    ind.vars = dofs_vars;
    ind.score = result_after_frags;
    //return ind.score;
    stats.increment("improved_after_de");
    stats.record("improved_amt_de",  std::abs( compare_antes - compare_result) );
    //result = result_after_frags;
    // core::pose::PoseOP check_pose = native->clone();
    // std::cout << "check pose dockindens" << std::endl;
    //ind.score = result;
    // std::cout << "init size " << size_pose << std::endl;
    // std::cout << "check_pose size " << check_pose->total_residue() << std::endl;
    // std::cout << "pose_with_frags size " << pose_with_frags->total_residue() << std::endl;
    // core::scoring::calpha_superimpose_pose( pose1, pose2 );
    // for ( int i = 1; i < inner_pose->total_residue(); ++i)  {
    // 	core::Size res_num_1 = core::Size(i);
    // 	core::PointPosition calpha1_pos  = inner_pose->residue(res_num_1).xyz("CA");
    // 	check_pose->residue(res_num_1).set_xyz(std::string("CA").c_str(), calpha1_pos);
    // }
    // double hcut_ = 0.18;
    // double pcut_ = 0.81;
    // double length_ = 38;
    // protocols::hybridization::DDomainParse ddom(pcut_,hcut_,length_);
    //utility::vector1<protocols::loops::Loops> domains = ddom.split( *inner_pose );
    // core::id::AtomID_Map< core::id::AtomID > atom_map;
    // core::pose::initialize_atomid_map( atom_map, *inner_pose, core::id::AtomID::BOGUS_ATOM_ID() );
    // protocols::hybridization::natom_aligned(*check_pose,*inner_pose,atom_map,1.0);

    // HybridMover->hybrid_protocol->align_by_domain(*check_pose, *inner_pose, domains);
    // core::Size atomno(0);
    // core::Vector x2;
    // for ( core::Size i=1; i<= 10; ++i ) {
    //   for ( core::Size j=1; j<= check_pose->residue_type(i).natoms(); ++j ) { // use residue_type to prevent internal coord update
    // 	++atomno;
    // 	x2 = inner_pose->xyz(core::id::AtomID(j, i) );
    // 	check_pose->set_xyz(core::id::AtomID( j,i), x2 );
    //   }
    // }

    // check_pose->residue(1);
    // check_pose = native->clone();
    // PoseDensityFragmentFunction::fill_pose(check_pose, ind, ss);
    // check_pose->residue(1);
    // check_pose->set_phi(1, 0);
    // check_pose->set_psi(check_pose->total_residue() - 1, 0);
    // if ( size_pose < check_pose->total_residue()) {
    //   core::pose::PoseOP recortado = core::pose::PoseOP(new core::pose::Pose(*check_pose,1, size_pose));
    //   check_pose = recortado->clone();
    // }

    // if (check_pose->is_centroid()) {
    // } else {
    //   protocols::moves::MoverOP tocen( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );
    //   tocen->apply(*check_pose);
    // }

    // for (int i = 1; i <= check_pose->total_residue(); i++) {
    //   check_pose->set_phi(i, pose_with_frags->phi(i));
    //   check_pose->set_psi(i, pose_with_frags->psi(i));
    //   check_pose->set_omega(i, pose_with_frags->omega(i));
    //   check_pose->set_secstruct(i,pose_with_frags->secstruct(i));
    // }
    // std::cout << "check_pose size " << check_pose->total_residue() << std::endl;
    // std::cout << "pose_with_frags size " << pose_with_frags->total_residue() << std::endl;

    // std::cout<< "sobreeescribe checkpose " << std::endl;
    // core::scoring::calpha_superimpose_pose(*check_pose, *pose_with_frags);
    // check_pose->residue(1);
    // dockindens->apply(*check_pose);
    // double result_3 = SCORE_ERROR_FIXED + (*scorefxn)(*check_pose);
    // std::cout << "RECALCULATE (ind) : " <<  (-1 * SCORE_ERROR_FIXED) + result_3  << std::endl;
    // check_pose->dump_pdb("./check_pose_superimpose.pdb");
    // pose_with_frags->dump_pdb("./frags_result_pose_Psuperimpose.pdb");
  }
  //  else {
  //   pose_to_ind(pose_backup, ind);
  //   ind.score = result;
  // }

  // std::cout << "FINAL IND SCORE " << ind.score << std::endl;
  result = ind.score;

  double end_score = (-1 * SCORE_ERROR_FIXED) +  ind.score;

  stats.record("avg_after_de", (-1 * SCORE_ERROR_FIXED) +  ind.score );

  stats.record("total_tries_de", stats.obtain("total_tries") );
  stats.erase("total_tries");

  stats.record("accepted_fragments_de", stats.obtain("accepted_fragments") );
  stats.erase("accepted_fragments");

  pose_ = native->clone();

  return result;
}

/*
    // for (int i = 1; i < check_pose->total_residue(); i++) {
    //   if ( std::abs(check_pose->phi(i) - pose_with_frags->phi(i)) > 0.1 ) {
    // 	std::cout << "different phi " << i << " vars " << check_pose->phi(i) << " != " << pose_with_frags->phi(i)
    // 		  << std::endl;
    //   }
    //   if ( std::abs(check_pose->psi(i) - pose_with_frags->psi(i)) > 0.1 ) {
    // 	std::cout << "different psi " << i << " vars " << check_pose->psi(i) << " != " << pose_with_frags->psi(i)
    // 		  << std::endl;
    //   }
    // }


 */
