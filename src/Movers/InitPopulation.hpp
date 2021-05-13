
#ifndef DIFFERENTIALEVOLUTION_INITPOPUL_HPP
#define DIFFERENTIALEVOLUTION_INITPOPUL_HPP

#include "../Algorithm/DE_types.hpp"
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
  virtual void disturb_individual(Individual &ind, int ind_size) {
  }
  
  virtual void apply(std::vector<Individual>& popul, int NP, int ind_size);

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
  boost::shared_ptr<protocols::simple_moves::CopyDofMover> copy_dofs;
  utility::vector1< std::pair< core::id::DOF_ID, core::Real > > vector_info_native;


  PosePopulation(FitFunctionPtr scfxn_ini , core::pose::PoseOP p) : InitPopulation(scfxn_ini) {
    pose_ = p->clone();
    dist_degrees = 5.0;
    std::map<core::Size, core::Size > res_map;
    for (int i = 1; i <= pose_->total_residue() - 1; i++) {
      res_map[i] = i;
    }

    copy_dofs = boost::shared_ptr<protocols::simple_moves::CopyDofMover>( new protocols::simple_moves::CopyDofMover( *pose_, res_map));
    //to_cen.apply(pose_);
    copy_dofs->apply( *pose_ );
    core::pose::copydofs::CopyDofsInfo dofs_info_container = copy_dofs->copy_dofs_info( *pose_ );
    vector_info_native =  dofs_info_container.dofs_info();
  }


  // PosePopulation(FitFunctionPtr scfxn_ini , core::pose::PoseOP p, std::string ss_) : InitPopulation(scfxn_ini) {
  //   pose_ = p->clone();
  //   dist_degrees = 5.0;
  //   ss = ss_;
  // }

  PosePopulation(FitFunctionPtr scfxn_ini , core::pose::PoseOP p, std::string ss_, double disturb_size) : InitPopulation(scfxn_ini) {
    pose_ = p->clone();
    dist_degrees = disturb_size;
    ss = ss_;
    std::map<core::Size, core::Size > res_map;
    for (int i = 1; i <= pose_->total_residue() - 1; i++) {
      res_map[i] = i;
    }

    copy_dofs = boost::shared_ptr<protocols::simple_moves::CopyDofMover>( new protocols::simple_moves::CopyDofMover( *pose_, res_map));
    //to_cen.apply(pose_);
    copy_dofs->apply( *pose_ );
    core::pose::copydofs::CopyDofsInfo dofs_info_container = copy_dofs->copy_dofs_info( *pose_ );
    vector_info_native =  dofs_info_container.dofs_info();
  }


  // PosePopulation(FitFunctionPtr scfxn_ini , core::pose::PoseOP p, double disturb_size) : InitPopulation(scfxn_ini) {
  //   pose_ = p->clone();
  //   dist_degrees = disturb_size;
  // }



  void apply(std::vector<Individual>& popul, int NP, int ind_size) ;
  virtual void disturb_individual(Individual &ind, int ind_size);
  void secstruct_individual(Individual &ind, int ind_size) ;
  void rama_individual(Individual &ind, int ind_size) ;
  double scale_rama( double old_value) ;
  double scale( double old_value) ;

public:
  core::pose::PoseOP pose_;
};
class PoseTotalRandomInitPopulation : public PosePopulation
{
public:

  PoseTotalRandomInitPopulation(FitFunctionPtr scfxn_ini , core::pose::PoseOP p, std::string ss_ ) : PosePopulation(scfxn_ini, p)  {
  }


  void disturb_individual(Individual &ind, int ind_size) ;
  void apply(std::vector<Individual>& popul, int NP, int ind_size) ;
};

class TwoStagesInitPopulation : public PosePopulation
{
public:
  double dist_degrees;
  std::string ss;
  core::pose::PoseOP my_pose_;
  FitFunctionPtr scorefxn;
  boost::shared_ptr<InitStagesMover> two_stages_mover;

  TwoStagesInitPopulation(FitFunctionPtr scfxn_ini , core::pose::PoseOP p, boost::shared_ptr<InitStagesMover> two_stages_mover_in, std::string ss_ ) : PosePopulation(scfxn_ini, p)  {
    my_pose_ = p;
    ss = ss_;
    scorefxn = scfxn_ini;
    two_stages_mover = two_stages_mover_in;
  }

  void disturb_individual(Individual &ind, int ind_size) ;

  void apply(std::vector<Individual>& popul, int NP, int ind_size) ;
};

#endif // DIFFERENTIALEVOLUTION_INITPOPUL_HPP




