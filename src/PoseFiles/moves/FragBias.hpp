
#ifndef FRAG_BIAS_HPP
#define FRAG_BIAS_HPP

#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <core/util/SwitchResidueTypeSet.hh>

//Density
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/elec/FA_ElecEnergy.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/lkball/LK_BallEnergy.hh>
#include <protocols/hybridization/util.hh>


//Cartesian sampler
#include <protocols/hybridization/CartesianSampler.hh>
#include <core/pose/selection.hh>
#include <protocols/hybridization/FragmentBiasAssigner.hh>

#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>

#include<core/fragment/FrameIterator.hh>

#include "PoseFunction.hpp"

class FragBiasCalculator
{
public:
  FragBiasCalculator() {
  }


  void apply(  core::pose::Pose &pose_ ) {

    protocols::hybridization::FragmentBiasAssigner frag_bias_assigner( pose_ );
   frag_bias_assigner.automode( pose_ , -0.5);

   core::fragment::FragSetOP frags = protocols::hybridization::create_fragment_set_no_ssbias(pose_, 5, 1, true);
     utility::vector1<core::fragment::FragSetOP> fragments_;
     fragments_.push_back(frags);

     utility::vector1<numeric::random::WeightedSampler> frag_bias_;
     frag_bias_ =  my_compute_frag_bias( pose_, fragments_ , frag_bias_assigner.getFragmentProbs()  );
  }


utility::vector1<numeric::random::WeightedSampler>
my_compute_frag_bias(
	core::pose::Pose &pose,
	utility::vector1<core::fragment::FragSetOP> fragment_sets,
	utility::vector1<core::Real> fragmentProbs_
){
        utility::vector1<numeric::random::WeightedSampler> frag_bias_tmp;
      	std::cout << "test this " << std::endl;
	// clean up the vector
	int wdw_to_freeze_ = 0;
	//utility::vector1<core::Real> fragmentProbs_;
	frag_bias_tmp.resize( fragment_sets.size(), 0.0 );

	// for each fragment size, smooth over the fragment window
	//   - handle mapping from frame->seqpos
	//   - don't allow any insertions that cross cuts
	for ( core::Size i_frag_set = 1; i_frag_set<=fragment_sets.size(); ++i_frag_set ) {
		utility::vector1< core::Real > frame_weights( pose.size(), 0.0 );

		for ( core::Size i_frame = 1; i_frame <= fragment_sets[i_frag_set]->nr_frames(); ++i_frame ) {
			core::fragment::ConstFrameIterator frame_it = fragment_sets[i_frag_set]->begin(); // first frame of the fragment library
			std::advance(frame_it, i_frame-1);  // point frame_it to the i_frame of the library
			core::Size seqpos_start = (*frame_it)->start();  // find starting and ending residue seqpos of the inserted fragment
			core::Size seqpos_end   = (*frame_it)->end();

			frame_weights[seqpos_start] = 0;
			for ( int i_pos = (int)seqpos_start; i_pos<=(int)seqpos_end; ++i_pos ) {
				frame_weights[seqpos_start] += fragmentProbs_[i_pos];
			}
			frame_weights[seqpos_start] /= (seqpos_end-seqpos_start+1);

			// don't allow insertions at cut points
			for ( int i_pos = (int)seqpos_start; i_pos<=(int)seqpos_end-1; ++i_pos ) {
				if ( pose.fold_tree().is_cutpoint(i_pos) ) {
					for ( int i_wdw = 0; i_wdw<wdw_to_freeze_; ++i_wdw ) {
						frame_weights[seqpos_start-i_wdw] = 0.0; // don't allow insertions at a frame where there are cut points downstream
						std::cout << "chainbreak at: " <<  i_pos << " seqpos_start: " << seqpos_start << " residue_to_freeze: " << seqpos_start-i_wdw << " wdw_to_freeze_: " << wdw_to_freeze_ << std::endl;
					}
				}
			} // assign
		} // each frame in a fragment set
		frag_bias_tmp[i_frag_set].weights(frame_weights);

		//////////////////////////////////////////////
		/////////////// for debug only ///////////////
		// report the probability for each frame in a fragment set
		for ( core::Size i_frame = 1; i_frame <= fragment_sets[i_frag_set]->nr_frames(); ++i_frame ) {
			core::fragment::ConstFrameIterator frame_it = fragment_sets[i_frag_set]->begin(); // first frame of the fragment library
			std::advance(frame_it, i_frame-1);  // point frame_it to the i_frame of the library
			core::Size seqpos_start = (*frame_it)->start();  // find starting and ending residue seqpos of the inserted fragment
			core::Size seqpos_end   = (*frame_it)->end();

			for ( int i_pos = (int)seqpos_start; i_pos<=(int)seqpos_end-1; ++i_pos ) { //
			  std::cout << "Prob_dens( " << i_frame << " " << seqpos_start << " ) = " << frame_weights[seqpos_start] << std::endl;
			}

			std::cout << "i_frame: " << i_frame << " seqpos_start " << seqpos_start <<  " seqpos_end " << seqpos_end  <<  " frame_weights: " <<  frame_weights[seqpos_start] << std::endl;

		} // frame /////////////// for debug only ///////////////
	} // each fragment set

	return frag_bias_tmp;
}
};



#endif // FRAG_BIAS_HPP
