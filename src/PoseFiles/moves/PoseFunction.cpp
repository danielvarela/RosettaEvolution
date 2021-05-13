
#include "PoseFunction.hpp"

#include "DE_types.hpp"
#include "FitFunction.hpp"

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>


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

#include "SecondStage.hpp"


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

