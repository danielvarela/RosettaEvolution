
#ifndef DIFFERENTIALEVOLUTION_INITPOPUL_HPP
#define DIFFERENTIALEVOLUTION_INITPOPUL_HPP

#include "DE_types.hpp"
#include "FitFunction.hpp"
#include "FragInsertionMover.hpp"
#include "PoseFunction.hpp"


#include <boost/shared_ptr.hpp>
#include <core/pose/Pose.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/RamaPrePro.hh>
#include <core/scoring/ScoringManager.hh>

class InitPopulation
{
public:
  InitPopulation(FitFunctionPtr scfxn_in) {
    scfxn = scfxn_in;
  }

  virtual void apply(std::vector<Individual>& popul, int NP, int ind_size) {

    for (int i = 0; i < NP; ++i) {
      Individual ind(ind_size);

      for (int j = 0; j < ind_size; ++j) {
	if (URAND() < 0.5) {
	  ind.vars[j] = URAND();
	} else {
	  ind.vars[j] = URAND() * -1;
	}
      }

      ind.vars[0] = 0;
      ind.vars[ind.vars.size() - 1] = 0;

      score(ind);

      popul.push_back(ind);
    }
  }

  double score(Individual& ind) {
     return scfxn->score(ind);
  }
protected:
  FitFunctionPtr scfxn;
};

typedef boost::shared_ptr<InitPopulation> InitPopulationPtr;

class PosePopulation : public InitPopulation
{
public:
  double dist_degrees;
  std::string ss;

  PosePopulation(FitFunctionPtr scfxn_ini , core::pose::PoseOP p) : InitPopulation(scfxn_ini) {
    pose_ = p->clone();
    dist_degrees = 150.0;
  }


  PosePopulation(FitFunctionPtr scfxn_ini , core::pose::PoseOP p, std::string ss_) : InitPopulation(scfxn_ini) {
    pose_ = p->clone();
    dist_degrees = 150.0;
    ss = ss_;
  }


  PosePopulation(FitFunctionPtr scfxn_ini , core::pose::PoseOP p, double disturb_size) : InitPopulation(scfxn_ini) {
    pose_ = p->clone();
    dist_degrees = disturb_size;
  }



  void apply(std::vector<Individual>& popul, int NP, int ind_size) {
    //    std::cout << "ind_size " <<  ind_size << std::endl;
    popul.resize(0);
    for (int i = 0; i < NP; ++i) {
      Individual ind(ind_size);
      disturb_individual(ind, ind_size);
      ind.vars[0] = 0;
      ind.vars[ind.vars.size() - 1] = 0;
      score(ind);
      popul.push_back(ind);
    }
  }


  virtual void disturb_individual(Individual &ind, int ind_size) {
    int nres = scfxn->start_res();
    dist_degrees = 180;
    for (int j = 0; j < ind_size - 2; j = j+2) {
      //taking into account the previous obtained angles
      if (URAND() > 0.5) {
	//ind.vars[j] = scale( std::min(180.0,  pose_->phi(nres) + (URAND() * dist_degrees) ) );
	ind.vars[j] =  pose_->phi(nres) + (URAND() * dist_degrees) ;

      } else {

	ind.vars[j] = pose_->phi(nres) - (URAND() * dist_degrees);
      }
        if (ind.vars[j] > 180.0) {
	  ind.vars[j] -= std::abs(ind.vars[j] - 180.0);
	}
	if (ind.vars[j] < -180.0) {
	  ind.vars[j] += std::abs(ind.vars[j] - 180.0);
	}
	ind.vars[j] = scale(ind.vars[j]);
      if (URAND() > 0.5) {
	ind.vars[j + 1] = pose_->psi(nres) + (URAND() * dist_degrees) ;
      } else {
	ind.vars[j + 1] = pose_->psi(nres) - (URAND() * dist_degrees);
      }
        if (ind.vars[j + 1] > 180.0) {
	  ind.vars[j + 1] -= std::abs(ind.vars[j + 1] - 180.0);
	}
	if (ind.vars[j + 1] < -180.0) {
	  ind.vars[j + 1] += std::abs(ind.vars[j + 1] - 180.0);
	}

	ind.vars[j+1] = scale(ind.vars[j+1]);

      nres++;
    }

    ind.omega.resize(0);
    for (int i = 1; i <= pose_->total_residue(); i++) {
      ind.omega.push_back(pose_->omega(i));
    }
    ind.ss = pose_->secstruct();
    ind.score = 1000;
  };

  void secstruct_individual(Individual &ind, int ind_size) {

    std::string ss_struct = pose_->secstruct();
    int nres = scfxn->start_res();
    for (int j = 0; j < ind_size - 2; j = j+2) {
      // secstruct phi and psi
      if (ss_struct[nres - 1] == 'H') {
	if (URAND() > 0.5) {
	  ind.vars[j] = scale(-64 + (URAND() * 7));
	} else {
	  ind.vars[j] = scale(-64 - (URAND() * 7));
	}

	if (URAND() > 0.5) {
	  ind.vars[j + 1] = scale(-41 + (URAND() * 7)  );
	} else {
	  ind.vars[j + 1] = scale(-41 - (URAND() * 7)  );
	}

      } else {

        if (URAND() < 0.5) {
	  ind.vars[j] = URAND();
	} else {
	  ind.vars[j] = URAND() * -1;
	}

      }


      nres++;
    }

  }

  void rama_individual(Individual &ind, int ind_size) {
    core::scoring::Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();

    std::string ss_struct = pose_->secstruct();

    double rama_phi, rama_psi;

    int nres = scfxn->start_res();
    for (int j = 0; j < ind_size - 2; j = j+2) {

      // rama phi and psi
      rama.random_phipsi_from_rama(pose_->aa(nres), rama_phi, rama_psi);
      ind.vars[j] = scale_rama(rama_phi);
      ind.vars[j + 1] = scale_rama(rama_psi);
      //	std::cout << "nres " << nres << " j " << j << " | " << rama_phi << " " << rama_psi << std::endl;
      nres++;
    }
  }

  double scale_rama( double old_value) {
    double old_min = 0.0;
    double old_max = 360.0;
    double new_min = -1*DE_LIMIT;
    double new_max = DE_LIMIT;

    return ( (old_value - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min;
  }

  double scale( double old_value) {
    double old_min = -180.0;
    double old_max = 180.0;
    double new_min = -1.0*DE_LIMIT;
    double new_max = DE_LIMIT;

    return ( (old_value - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min;
  }

public:
  core::pose::PoseOP pose_;
};
class PoseTotalRandomInitPopulation : public PosePopulation
{
public:

  PoseTotalRandomInitPopulation(FitFunctionPtr scfxn_ini , core::pose::PoseOP p, std::string ss_ ) : PosePopulation(scfxn_ini, p)  {
  }


  void disturb_individual(Individual &ind, int ind_size) {
     for (int j = 0; j < ind_size; ++j) {
	if (URAND() < 0.5) {
	  ind.vars[j] = URAND();
	} else {
	  ind.vars[j] = URAND() * -1;
	}
      }

      ind.vars[0] = 0;
      ind.vars[ind.vars.size() - 1] = 0;
      ind.omega.resize(0);
      for (int i = 1; i <= pose_->total_residue(); i++) {
	ind.omega.push_back(pose_->omega(i));
      }
      ind.ss = pose_->secstruct();
      ind.score = 1000;
  }
  void apply(std::vector<Individual>& popul, int NP, int ind_size) {

    boost::shared_ptr<PoseFragmentFunction> pfunc = boost::dynamic_pointer_cast<PoseFragmentFunction >(scfxn);
    for (int i = 0; i < NP; ++i) {
      Individual ind(ind_size);
      disturb_individual(ind, ind_size);
      pfunc->fill_pose(pose_, ind, ss);
      pfunc->pose_to_ind(pose_, ind );
      score(ind);
      popul.push_back(ind);
    }

  }
};

class TwoStagesInitPopulation : public PosePopulation
{
public:
  double dist_degrees;
  std::string ss;
  core::pose::PoseOP pose_;
  FitFunctionPtr scorefxn;
  boost::shared_ptr<InitStagesMover> two_stages_mover;

  TwoStagesInitPopulation(FitFunctionPtr scfxn_ini , core::pose::PoseOP p, boost::shared_ptr<InitStagesMover> two_stages_mover_in, std::string ss_ ) : PosePopulation(scfxn_ini, p)  {
    pose_ = p;
    ss = ss_;
    scorefxn = scfxn_ini;
    two_stages_mover = two_stages_mover_in;
  }


  void apply(std::vector<Individual>& popul, int NP, int ind_size) {

    boost::shared_ptr<PoseScoreFunction> pfunc = boost::dynamic_pointer_cast<PoseScoreFunction >(scorefxn);
    for (int i = 0; i < NP; ++i) {
      Individual ind(ind_size);
      disturb_individual(ind, ind_size);
      pfunc->fill_pose(pose_, ind, ss);
      two_stages_mover->apply(*pose_);
      pfunc->pose_to_ind(pose_, ind );
      ind.vars[0] = 0;
      ind.vars[ind.vars.size() - 1] = 0;
      score(ind);
      popul.push_back(ind);
    }

  }
};

#endif // DIFFERENTIALEVOLUTION_INITPOPUL_HPP




