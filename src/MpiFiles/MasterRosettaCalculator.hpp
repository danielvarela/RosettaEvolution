#ifndef MASTERROSETTACALCULATOR_H
#define MASTERROSETTACALCULATOR_H

#if(MPI_ENABLED)


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

namespace mpi = boost::mpi;


#include "../Movers/CalculateRmsdDistancePopulation.hpp"
#include <protocols/relax/FastRelax.hh>

class MasterRosettaCalculator {
public:
  mpi::communicator world;
  boost::shared_ptr<FitFunction> scfxn_ind_with_frags;
  boost::shared_ptr<PoseScoreFunction> simple_ffxn;
  boost::shared_ptr<PoseDensityFunction> density_func;
  std::string stage_name;
  explicit MasterRosettaCalculator(boost::shared_ptr<FitFunction> scfxn_ind_with_frags_input);

  protocols::electron_density::SetupForDensityScoringMoverOP dockindens;
  boost::shared_ptr<protocols::simple_moves::SwitchResidueTypeSetMover> to_all;
  boost::shared_ptr<protocols::simple_moves::SwitchResidueTypeSetMover> to_cen;
  core::scoring::ScoreFunctionOP fa_score;
  boost::shared_ptr<protocols::relax::FastRelax> my_relax;

  std::vector<std::vector<IndMPI> >
  SplitVector(const std::vector<Individual>& vec, size_t n);


  std::vector<Individual>
  run(std::vector<Individual> num_popul, int mode = 2);

  std::vector<Individual>
  dump(std::vector<Individual> num_popul);

  virtual std::vector<IndMPI>
  run(std::vector<Individual> popul, std::vector<Individual> trial_popul, int mode = 2) {
    std::cout << "error, wrong initialization of MasterRosetttaCalculator" << std::endl;
    exit(1);
  }


};

class MasterEvaluateAndNearestCalculator : public MasterRosettaCalculator {
public:
  CalculateDistancePopulationPtr calculate_distances_popul;

  explicit MasterEvaluateAndNearestCalculator(boost::shared_ptr<FitFunction> scfxn_ind_with_frags_input,   CalculateDistancePopulationPtr  calculate_distances_popul_in) : MasterRosettaCalculator(scfxn_ind_with_frags_input) {
    calculate_distances_popul = calculate_distances_popul_in;
  }

  std::vector<IndMPI>
  run(std::vector<Individual> popul, std::vector<Individual> trial_popul, int mode = 2) override;

};


class MasterScatterGather : public MasterRosettaCalculator {
public:
  CalculateDistancePopulationPtr calculate_distances_popul;

  explicit MasterScatterGather(boost::shared_ptr<FitFunction> scfxn_ind_with_frags_input,   CalculateDistancePopulationPtr  calculate_distances_popul_in) : MasterRosettaCalculator(scfxn_ind_with_frags_input) {
    calculate_distances_popul = calculate_distances_popul_in;
  }

  std::vector<IndMPI>
  run(std::vector<Individual> popul, std::vector<Individual> trial_popul, int mode = 2) override;

};



#endif

#endif /* MASTERROSETTACALCULATOR_H */
