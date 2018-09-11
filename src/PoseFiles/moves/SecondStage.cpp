#include "SecondStage.hpp"

SecondStage::SecondStage() {
  to_fa = protocols::simple_moves::SwitchResidueTypeSetMoverOP(new protocols::simple_moves::SwitchResidueTypeSetMover("fa_standard"));
  to_cen =  protocols::simple_moves::SwitchResidueTypeSetMoverOP(new protocols::simple_moves::SwitchResidueTypeSetMover("centroid"));

  // 	// set up minimizer
  rbminimizer = core::optimization::AtomTreeMinimizerOP( new core::optimization::AtomTreeMinimizer  ) ;
  options_rb = core::optimization::MinimizerOptionsOP( new  core::optimization::MinimizerOptions("lbfgs_armijo_nonmonotone", 0.01, true, false, false) );
  options_rb->max_iter(20);
  mm_rb = core::kinematics::MoveMapOP( new core::kinematics::MoveMap);
  mm_rb->set_bb ( true );
  mm_rb->set_chi ( true );
  mm_rb->set_jump ( true );

  densonly = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() );
  densonly->set_weight( core::scoring::elec_dens_fast, 5.0 );

  //     	pack_mover = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover );
}

void
SecondStage::apply(core::pose::PoseOP p) {
  // std::cout << "to cen " << std::endl;
  to_cen->apply( *p );
  std::cout << "end to cen " << std::endl;
  std::cout << "score " << std::endl;
  (*densonly)(*p);
  rbminimizer->run( *p, *mm_rb, *densonly, *options_rb );
  to_fa->apply( *p );
  //    pack_mover->apply( *p );

}
