
#ifndef WORKERPROCESS_H
#define WORKERPROCESS_H



#include <iostream>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include "mpi.h"
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <vector>
#include <map>
#include <string>
#include "../Controller/DE_Operator.hpp"
#include "UtilsMPI.hpp"
#include <protocols/relax/FastRelax.hh>

namespace mpi = boost::mpi;

class WorkerProcess
{
public:
  boost::shared_ptr<FitFunction> scfxn_ind_with_frags_1,  scfxn_ind_simple_1;
  boost::shared_ptr<FitFunction> scfxn_ind_with_frags_2,  scfxn_ind_simple_2;
  boost::shared_ptr<FitFunction> scfxn_ind_with_frags_3,  scfxn_ind_simple_3;
  boost::shared_ptr<FitFunction> scfxn_ind_with_frags_complete,  scfxn_ind_simple_complete;
  boost::shared_ptr<FitFunction> scfxn_density;


  protocols::electron_density::SetupForDensityScoringMoverOP dockindens;
  boost::shared_ptr<protocols::simple_moves::SwitchResidueTypeSetMover> to_all;
  boost::shared_ptr<protocols::simple_moves::SwitchResidueTypeSetMover> to_cen;
  core::scoring::ScoreFunctionOP fa_score;
  boost::shared_ptr<protocols::relax::FastRelax> my_relax;

  mpi::communicator world;

  class ScoreStrategy {
  public:
    boost::shared_ptr<FitFunction> scfxn_ind_with_frags_1,  scfxn_ind_simple_1;
    boost::shared_ptr<FitFunction> scfxn_ind_with_frags_2,  scfxn_ind_simple_2;
    boost::shared_ptr<FitFunction> scfxn_ind_with_frags_3,  scfxn_ind_simple_3;
    boost::shared_ptr<FitFunction> scfxn_ind_with_frags_complete,  scfxn_ind_simple_complete;
    boost::shared_ptr<FitFunction> scfxn_density;
    boost::shared_ptr<PoseDensityFunction> scfxn_density_simple;
    int mode;
    std::string stage_name;
    ScoreStrategy(   boost::shared_ptr<FitFunction> scfxn_ind_with_frags_1_in,  boost::shared_ptr<FitFunction> scfxn_ind_simple_1_in,
		     boost::shared_ptr<FitFunction> scfxn_ind_with_frags_2_in, boost::shared_ptr<FitFunction>  scfxn_ind_simple_2_in,
		     boost::shared_ptr<FitFunction> scfxn_ind_with_frags_3_in, boost::shared_ptr<FitFunction>  scfxn_ind_simple_3_in,
		     boost::shared_ptr<FitFunction> scfxn_ind_with_frags_complete_in, boost::shared_ptr<FitFunction>  scfxn_ind_simple_complete_in,
		     int mode_in, std::string stage_name_in
		     );
    void apply(Individual& ind);
  };

  explicit WorkerProcess(boost::shared_ptr<DE_Operator> app_operator);
  void run();

};

class WorkerProcessEvaluateAndNearest : public WorkerProcess
{
public:
  CalculateDistancePopulationPtr calculate_distances_popul;

  explicit WorkerProcessEvaluateAndNearest(boost::shared_ptr<DE_Operator> app_operator);


  void run();
};


class WorkerProcessScatterGather : public WorkerProcess
{
public:
  CalculateDistancePopulationPtr calculate_distances_popul;
  explicit WorkerProcessScatterGather(boost::shared_ptr<DE_Operator> app_operator);

  void run();
};



#endif /* WORKERPROCESS_H */
