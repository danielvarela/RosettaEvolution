#ifndef DIFFERENTIALEVOLUTION_POSEFUNCTION_HPP
#define DIFFERENTIALEVOLUTION_POSEFUNCTION_HPP

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/copydofs/util.hh>
#include <core/pose/copydofs/CopyDofs.hh>
#include <core/pose/copydofs/CopyDofsInfo.hh>
#include <protocols/simple_moves/CopyDofMover.hh>
#include <core/pose/copydofs/util.hh>


//Density
#include <protocols/electron_density/SetupForDensityScoringMover.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/fragment/Frame.hh>

#include "../Algorithm/DE_types.hpp"
#include "FitFunction.hpp"
#include "FragInsertionMover.hpp"


class PoseFunction : public FitFunction
{
public:
  core::pose::PoseOP native;
  PoseFunction() {}
  PoseFunction(core::pose::PoseOP p) : FitFunction() {
    pose_ = p->clone();
    native = p->clone();
    scorefxn = core::scoring::get_score_function();
    //scorefxn->set_weight( core::scoring::elec_dens_window, 1.0 );
    D_ = (pose_->size()) * 2;
    //D_ = pose_->size() * 2;
  }

  PoseFunction(core::pose::PoseOP p , core::scoring::ScoreFunctionOP sfxn, std::string ss_in) : FitFunction() {
    pose_ = p->clone();
    scorefxn = sfxn;
    D_ = (pose_->size()) * 2;
    //    D_ = pose_->size() * 2;
    //    std::cout << "D_ constructor " << D_ << std::endl;
    ss = ss_in;
    native = p->clone();

  }

  virtual double score(Individual& ind) override {
    core::pose::PoseOP inner_pose = native->clone();
    fill_pose(inner_pose, ind, ss);
    double result = (*scorefxn)(*inner_pose);
    ind.score = result;
    return result;
  }

  void fill_pose(core::pose::PoseOP& p, const Individual& ind, std::string ss_) {
    Individual ind_converted = convert(ind);
    int nres = start_res();
    int size = D_;
    p = native->clone();


    // std::cout << "fill_pose " << std::endl;
    for (int i = 0 ; i < size && (nres <= static_cast<int>(p->size())); ++i) {
      //  std::cout << nres << " : " << p->sequence()[nres - 1] << std::endl;
      if (is_pair(i)) {
	p->set_phi(nres, ind_converted.vars[i]);
      } else {
        p->set_psi(nres, ind_converted.vars[i]);
        p->set_secstruct(nres, ind_converted.ss[nres - 1]);
	nres++;
      }
    }

    for (int i = 1; i < p->total_residue(); i++) {
      p->set_omega(i, ind_converted.omega[i - 1]);
    }


  }

  int start_res() {
    return 1;
  }

  static bool is_pair(int i) {
      if ( i % 2 == 0) {
	return true;
	  } else {
        return false;
      }
  }


  int D() {
    return D_;
  }

  double lim() {
    return 180.0;
  }

  std::string name() {
    return name_;
  }

public:
  int D_;
  core::pose::PoseOP pose_;
  core::scoring::ScoreFunctionOP scorefxn;
  std::string ss;

};

class PoseScoreFunction : public PoseFunction
{
public:
  FragInsertionMoverPtr frag_mover;


  PoseScoreFunction(){}
  PoseScoreFunction(core::pose::PoseOP p, core::scoring::ScoreFunctionOP sfxn, std::string ss_in, FragInsertionMoverPtr frag_mover_) : PoseFunction(p, sfxn, ss_in) {
    frag_mover = frag_mover_;
  }

  double scale( double old_value) {
    double old_min = -180.0;
    double old_max = 180.0;
    double new_min = -1 * DE_LIMIT;
    double new_max = DE_LIMIT;

    return ( (old_value - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min;
  }

  std::string print_stats() override {
    int total_improves = 0;
    total_improves = static_cast<int>(stats.obtain("improved_after_de"));

    std::string output_string = "[DE_IMPROVED] " + std::to_string(total_improves);
    if (total_improves > 0) {
      output_string += " amt_improves " + std::to_string( stats.obtain("improved_amt_de") );
    } else {
      output_string += " amt_improves 0";
    }
    if (stats.obtain("avg_before_de") != 0) {
      output_string += " avg_before_de " + std::to_string(stats.obtain("avg_before_de") / stats.obtain("de_individuals") );
    }
    if (stats.obtain("avg_after_de") != 0) {
      output_string += " avg_after_de " + std::to_string(stats.obtain("avg_after_de") / stats.obtain("de_individuals") );
    }


    if (stats.obtain("total_tries_de") > 0)
      output_string += " total_tries " + std::to_string(stats.obtain("total_tries_de") );

    if (stats.obtain("accepted_fragments_de") > 0)
      output_string += " accepted_fragments " + std::to_string(stats.obtain("accepted_fragments_de"));

    return output_string;
  }

  void fill_pose(core::pose::PoseOP& p, const Individual& ind, std::string ss_) {
    Individual ind_converted = convert(ind);
    int nres = start_res();
    int size = D_;
    p = native->clone();


    for (int i = 0 ; i < size && nres <= p->total_residue() - 1; i=i+2) {
	p->set_phi(nres, ind_converted.vars[i]);
        p->set_psi(nres, ind_converted.vars[i+1]);
	nres++;
      }

    for (int i = 1; i < p->total_residue() - 1; i++) {
      p->set_secstruct(i, ind_converted.ss[i-1]);
    }


    for (int i = 1; i < p->total_residue() - 1; i++) {
      p->set_omega(i, ind_converted.omega[i - 1]);
    }

  }


  void pose_to_ind(core::pose::PoseOP inner_pose_, Individual& ind) {
    int ind_size = D_;
    int nres = 1;
    for (int j = 0; nres < inner_pose_->total_residue(); j = j+2) {
      ind.vars[j] = scale(inner_pose_->phi(nres));
      ind.vars[j + 1] = scale(inner_pose_->psi(nres));
      ind.omega[nres - 1] = inner_pose_->omega(nres);
      nres++;
    }
    ind.ss = inner_pose_->secstruct();
  }

  virtual double score(Individual& ind) override {
    core::pose::PoseOP inner_pose = native->clone();
    fill_pose(inner_pose, ind, ss);
    double result = SCORE_ERROR_FIXED + (*scorefxn)(*inner_pose);
    ind.score = result;
    return result;
   }

  double apply_fragment_stage(Individual& ind) {
    return ind.score;
  }

  virtual double run_frag_mover(core::pose::Pose& pose) {
    frag_mover->apply(pose);
    return SCORE_ERROR_FIXED +  (*scorefxn)(pose);
  }

};

class PoseDensityFunction : public PoseScoreFunction
{
public:

  core::pose::PoseOP native;
  protocols::electron_density::SetupForDensityScoringMoverOP dockindens;
  boost::shared_ptr<protocols::simple_moves::SwitchResidueTypeSetMover> to_all;
  boost::shared_ptr<protocols::simple_moves::CopyDofMover> copy_dofs;
  utility::vector1< std::pair< core::id::DOF_ID, core::Real > > vector_info_native;

  PoseDensityFunction(core::pose::PoseOP p, core::scoring::ScoreFunctionOP sfxn, std::string ss_in, FragInsertionMoverPtr frag_mover_) : PoseScoreFunction(p, sfxn, ss_in, frag_mover_) {
    native =  p->clone();
    dockindens = protocols::electron_density::SetupForDensityScoringMoverOP( new protocols::electron_density::SetupForDensityScoringMover );
    // scorefxn = core::scoring::ScoreFunctionFactory::create_score_function("beta_cart");
    scorefxn->set_weight( core::scoring::elec_dens_fast, 30.0 );
    to_all = boost::shared_ptr<protocols::simple_moves::SwitchResidueTypeSetMover>( new protocols::simple_moves::SwitchResidueTypeSetMover("fa_standard"));
    std::map<core::Size, core::Size > res_map;
    dockindens->apply(*native);
    for (int i = 1; i <= native->total_residue() - 1; i++) {
      res_map[i] = i;
    }
    copy_dofs = boost::shared_ptr<protocols::simple_moves::CopyDofMover>( new protocols::simple_moves::CopyDofMover( *native, res_map));
    //to_cen.apply(pose_);
    copy_dofs->apply( *native );
    core::pose::copydofs::CopyDofsInfo dofs_info_container = copy_dofs->copy_dofs_info( *native );
    vector_info_native =  dofs_info_container.dofs_info();
    double dofs_size = vector_info_native.size();
    D_ = dofs_size;
  }

  virtual double score(Individual& ind) override;

  virtual double run_frag_mover(core::pose::Pose& pose);

  void fill_pose(core::pose::PoseOP& p, const Individual& ind, std::string ss_) override;
};


class PoseFragmentFunction : public PoseScoreFunction {
public:
  PoseFragmentFunction() : PoseScoreFunction() {}
  PoseFragmentFunction(core::pose::PoseOP p, core::scoring::ScoreFunctionOP sfxn, std::string ss_in, FragInsertionMoverPtr frag_mover_) : PoseScoreFunction(p, sfxn, ss_in, frag_mover_) {
  }

  double scale( double old_value) {
    double old_min = -180.0;
    double old_max = 180.0;
    double new_min = -1 * DE_LIMIT;
    double new_max = DE_LIMIT;

    return ( (old_value - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min;
  }

  std::string print_stats() override {
    int total_improves = 0;
    total_improves = static_cast<int>(stats.obtain("improved_after_de"));

    std::string output_string = "[DE_IMPROVED] " + std::to_string(total_improves);
    if (total_improves > 0) {
      output_string += " amt_improves " + std::to_string( stats.obtain("improved_amt_de") );
    } else {
      output_string += " amt_improves 0";
    }

    if (stats.obtain("avg_before_de") != 0) {
      output_string += " avg_before_de " + std::to_string(stats.obtain("avg_before_de") / stats.obtain("de_individuals") );
    }
    if (stats.obtain("avg_after_de") != 0) {
      output_string += " avg_after_de " + std::to_string(stats.obtain("avg_after_de") / stats.obtain("de_individuals") );
    }


    if (stats.obtain("total_tries_de") > 0.001) {
      output_string += " total_tries " + std::to_string(stats.obtain("total_tries_de") );
    }

    if (stats.obtain("accepted_fragments_de") > 0.001) {
    output_string += " accepted_fragments " + std::to_string(stats.obtain("accepted_fragments_de"));
    }
    return output_string;
  }

  double score(Individual& ind) override {
    return apply_fragment_stage(ind);
  }

  virtual double apply_fragment_stage(Individual& ind);

  virtual double run_frag_mover(core::pose::Pose& pose);

};

class PoseDensityFragmentFunction : public PoseFragmentFunction
{
public:
  core::pose::PoseOP native;
  protocols::electron_density::SetupForDensityScoringMoverOP dockindens;
  boost::shared_ptr<protocols::simple_moves::SwitchResidueTypeSetMover> to_all;
  protocols::simple_moves::SwitchResidueTypeSetMoverOP to_cen;
  boost::shared_ptr<protocols::simple_moves::CopyDofMover> copy_dofs;
  utility::vector1< std::pair< core::id::DOF_ID, core::Real > > vector_info_native;
  core::scoring::ScoreFunctionOP scorefxn_fa;

  PoseDensityFragmentFunction() : PoseFragmentFunction() {}
  PoseDensityFragmentFunction(core::pose::PoseOP p, core::scoring::ScoreFunctionOP sfxn, std::string ss_in, FragInsertionMoverPtr frag_mover_) :  PoseFragmentFunction(p, sfxn, ss_in, frag_mover_) {
  dockindens = protocols::electron_density::SetupForDensityScoringMoverOP( new protocols::electron_density::SetupForDensityScoringMover );
  scorefxn_fa = core::scoring::ScoreFunctionFactory::create_score_function("beta_cart");
  scorefxn->set_weight( core::scoring::elec_dens_fast, 30.0 );
  scorefxn_fa->set_weight( core::scoring::elec_dens_fast, 30.0 );
  to_all = boost::shared_ptr<protocols::simple_moves::SwitchResidueTypeSetMover>( new protocols::simple_moves::SwitchResidueTypeSetMover("fa_standard"));
  //to_cen.apply(pose_);
  to_cen =  protocols::simple_moves::SwitchResidueTypeSetMoverOP(new protocols::simple_moves::SwitchResidueTypeSetMover("centroid"));
    std::map<core::Size, core::Size > res_map;
  native = p->clone();
    dockindens->apply(*native);
    for (int i = 1; i <= native->total_residue() - 1; i++) {
      res_map[i] = i;
    }

  copy_dofs = boost::shared_ptr<protocols::simple_moves::CopyDofMover>( new protocols::simple_moves::CopyDofMover( *native, res_map));
  copy_dofs->apply( *native );

  //to_cen.apply(pose_);
  core::pose::copydofs::CopyDofsInfo dofs_info_container = copy_dofs->copy_dofs_info( *native );
  vector_info_native =  dofs_info_container.dofs_info();
  double dofs_size = vector_info_native.size();
  D_ = dofs_size;

  }
  void
  apply_relax_stage(core::pose::PoseOP & inner_pose);

  double apply_fragment_stage(Individual& ind) override;


  double run_frag_mover(core::pose::Pose& pose) {
    frag_mover->apply(pose);
    dockindens->apply(pose);
    double result = SCORE_ERROR_FIXED +  (*scorefxn)(pose);
    return result;
  }

  void fill_pose(core::pose::PoseOP& p, const Individual& ind, std::string ss_) override;
};

class PoseToyFunction : public PoseFragmentFunction
{
public:

  PoseToyFunction(core::pose::PoseOP p, core::scoring::ScoreFunctionOP sfxn, std::string ss_in, FragInsertionMoverPtr frag_mover_) : PoseFragmentFunction(p, sfxn, ss_in, frag_mover_) {
  }
  double score(Individual& ind) override {
    int size = D_;

    double acc = 0.0;
    for (int i = 0 ; i < size; ++i) {
      acc += std::pow(ind.vars[i],2);
    }

    acc = SCORE_ERROR_FIXED + acc;
    ind.score = acc;
    return acc;
  }

};

struct PoseRange
{
  int i, j;
};

class PoseRangeFunction : public PoseFunction
{
public:
  PoseRangeFunction(core::pose::PoseOP p, int i, int j, core::scoring::ScoreFunctionOP sfxn, std::string ss_in) : PoseFunction(p, sfxn, ss_in) {
    range.i = i;
    range.j = j;
    D_ = (j - i) * 2;

  }

  void fill_pose(core::pose::PoseOP& p, const Individual& ind, std::string ss_) {
    Individual ind_converted = convert(ind);
    int nres = range.i;
    int size = D_;
    p = native->clone();


    for (int i = 0 ; i < size  && (nres < static_cast<int>(p->size())); ++i) {
      if (is_pair(i)) {
	p->set_phi(nres, ind_converted.vars[i]);
      } else {
        p->set_psi(nres, ind_converted.vars[i]);
        p->set_secstruct(nres, ind_converted.ss[nres - 1]);
	nres++;
      }
    }
    for (int i = 1; i < p->total_residue(); i++) {
      p->set_omega(i, ind_converted.omega[i - 1]);
    }

  }

  int start_res() {
    return range.i;
  }

  std::string name() {
    return std::string("Pose Range Function");
  }

protected:
  PoseRange range;
};

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pack/pack_rotamers.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <core/conformation/util.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/conformation/util.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/SingleResidueFragData.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>

#include <core/fragment/BBTorsionAndAnglesSRFD.hh>

#include "SecondStage.hpp"

class PosePartialFunction : public PoseFunction
{
public:
  PosePartialFunction(core::pose::PoseOP p, int i, int j, core::scoring::ScoreFunctionOP sfxn, std::string ss_in) : PoseFunction(p, sfxn, ss_in) {
    range.i = i;
    range.j = j;
    D_ = (j - i) * 2;

    ppose = core::pose::PoseOP(new core::pose::Pose(*p, i, j));

    sstage = boost::shared_ptr<SecondStage >( new SecondStage() );
  std::cout << "start density align" << std::endl;


   dockindens = protocols::electron_density::SetupForDensityScoringMoverOP( new protocols::electron_density::SetupForDensityScoringMover );


  dockindens->apply( *ppose );
std::cout << "start score at res " << i << " to " << j << " : " <<  (*scorefxn)(*ppose) << std::endl;
}

  double score(Individual& ind) {
    Individual ind_converted = convert(ind);
    fill_pose(ppose, ind_converted, ss);
    dockindens->apply( *ppose );

    // sstage->apply(ppose);

    double result = (*scorefxn)(*ppose);

    ind.score = result;
    return result;
  }

  void fill_pose(core::pose::PoseOP& p, const Individual& ind, std::string ss_) {
    int nres = 1;
    if (static_cast<int>(p->size()) > D_) {
      nres = range.i;

      std::cout << "fill best pose from " << range.i << std::endl;
    }

    int size = D_;

    for (int i = 0 ; i < size; ++i) {
     if (is_pair(i) && (nres <= static_cast<int>(p->size()))) {
	p->set_phi(nres, ind.vars[i]);
      } else {
        p->set_psi(nres, ind.vars[i]);
        p->set_secstruct(nres, ind.ss[nres - 1]);
	nres++;
      }
    }

    for (int i = 1; i < p->total_residue(); i++) {
      p->set_omega(i, pose_->omega(i));
    }

  }

  int start_res() {
    return range.i;
  }

  std::string name() {
    return std::string("Pose Partial Function");
  }

protected:
  PoseRange range;
  core::pose::PoseOP ppose;
  boost::shared_ptr<SecondStage > sstage;
protocols::electron_density::SetupForDensityScoringMoverOP dockindens;
};

class PoseImprovedFunction : public PoseRangeFunction
{
public:

  PoseImprovedFunction(core::pose::PoseOP p, int i, int j, core::scoring::ScoreFunctionOP sfxn, std::string ss_in) : PoseRangeFunction(p, i, j, sfxn, ss_in) {

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


using namespace core::pack::task;
using core::pack::task::operation::TaskOperationCOP;


    // 	pack_mover = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover );
    // // 	// pack_mover->task_factory( main_task_factory );
    // 	// pack_mover->score_function( nonsymm_fa_scorefxn );

  }


  double score(Individual& ind);

    std::string name() {
    return std::string("Pose Improved Function");
  }

protected:
  protocols::simple_moves::SwitchResidueTypeSetMoverOP to_fa;
  protocols::simple_moves::SwitchResidueTypeSetMoverOP to_cen;
  core::optimization::AtomTreeMinimizerOP rbminimizer;
  core::optimization::MinimizerOptionsOP  options_rb ;
  core::kinematics::MoveMapOP mm_rb;
  core::scoring::ScoreFunctionOP densonly; 
};


#endif // DIFFERENTIALEVOLUTION_POSEFUNCTION_HPP
