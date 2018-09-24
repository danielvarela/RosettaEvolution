
#ifndef CALCULATERMSDDISTANCEPOPULATION_H
#define CALCULATERMSDDISTANCEPOPULATION_H

#include "DE_types.hpp"

#include <map>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include "FitFunction.hpp"
#include "PoseFunction.hpp"

class InterDistancesGraph
{
public:
  std::vector<double> bin_values;
  std::map<double, int> occurrence_map;
  InterDistancesGraph();

  void add_occurrence(double value);
  std::string print_graph();
};

class NeighStruct
{
public:
  NeighStruct() : index(0), distance(0) {}
  NeighStruct(int i, double d) : index(i), distance(d) {}

  int index;
  double distance;
};

class CalculateDistancePopulation
{
public:
  core::pose::PoseOP pose_ind, pose_other;
  std::vector<core::pose::PoseOP> popul_pdb;
  FitFunctionPtr scorefxn;
  std::string ss;
  core::pose::PoseOP native_;
  boost::shared_ptr<PoseFunction> pfunc;
  double fit_radius = 3.0;

  class DistancesResult
  {
  public:
    DistancesResult() {}

    std::vector<double> shared_fitness, rmsd_to_native, distances_of_population;
    std::vector<std::vector<NeighStruct> > neigh_per_ind;
  };


  CalculateDistancePopulation() {
  }

  CalculateDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius);


  DistancesResult
  run(const std::vector<Individual>& popul);

  virtual void
  apply_to_population(const std::vector<Individual>& popul, std::vector<double>& shared_fitness, std::vector<double>& rmsd_to_native,std::vector<double>& distances_of_population, std::vector<std::vector<NeighStruct> >& neigh_per_ind  );

  double distance_of_individual(core::pose::PoseOP pose_ind,const std::vector<core::pose::PoseOP>& popul_pdb, double& shared_acc, std::vector<NeighStruct>& neighs, std::vector<double>& ind_distances);

  void print_distance_population( const std::vector<Individual>& popul ) ;

  void
  calculate_rmsd(const std::vector<Individual>& popul,    std::vector<double>& rmsd_to_native );

  int
  find_nearest(Individual target, const std::vector<Individual>& popul);

  std::vector<int>
  find_nearest_parent(int ind_idx, Individual target, const std::vector<Individual>& popul);

  double distance_between_inds(Individual ind_left, Individual ind_right);

  virtual double
  current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2) = 0;

  std::vector<double>
  run_rmsd_to_native(const std::vector<Individual>& popul);

  void build_pdb_population(const std::vector<Individual>& population, std::vector<core::pose::PoseOP>& popul_pdb);

 };


class CalculateRmsdDistancePopulation : public CalculateDistancePopulation
{
public:
  CalculateRmsdDistancePopulation() : CalculateDistancePopulation() {

  }

  CalculateRmsdDistancePopulation(const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius) : CalculateDistancePopulation(p, sfxn, ss_in, pscore, radius) {}

  double
  current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2);
};

class CalculateNativeDiffDistancePopulation : public CalculateDistancePopulation
{
public:
  CalculateNativeDiffDistancePopulation() : CalculateDistancePopulation() {
  }

  CalculateNativeDiffDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius) ;

  double
  current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2);
};


class CalculateEuclideanDistancePopulation : public CalculateDistancePopulation
{
public:
  CalculateEuclideanDistancePopulation() : CalculateDistancePopulation() {
  }

  CalculateEuclideanDistancePopulation(FitFunctionPtr fitfunc);

  CalculateEuclideanDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius) ;

  double
  current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2);

  double
  euclidean_two_vectors(const std::vector<double> & vars_1, const std::vector<double>& vars_2);

};

class CalculateEuclideanLoopDistancePopulation : public CalculateEuclideanDistancePopulation
{
public:
  CalculateEuclideanLoopDistancePopulation() : CalculateEuclideanDistancePopulation() {
  }

  CalculateEuclideanLoopDistancePopulation(FitFunctionPtr fitfunc);

  CalculateEuclideanLoopDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius) ;

  double
  current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2);
};


class CalculateEuclideanDiffAbsDistancePopulation : public CalculateEuclideanDistancePopulation
{
public:
  CalculateEuclideanDiffAbsDistancePopulation() : CalculateEuclideanDistancePopulation() {
  }

  CalculateEuclideanDiffAbsDistancePopulation(FitFunctionPtr fitfunc);

  CalculateEuclideanDiffAbsDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius) ;

  double
  current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2);

};



class CalculateEuclideanMarioDistancePopulation : public CalculateEuclideanDistancePopulation
{
public:
  std::vector<int> selected_residues_for_rmsd;

  CalculateEuclideanMarioDistancePopulation() : CalculateEuclideanDistancePopulation() {
  }

  CalculateEuclideanMarioDistancePopulation(FitFunctionPtr fitfunc);

  CalculateEuclideanMarioDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius) ;

  double
  current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2);

};

class CalculateEuclideanMarioPartialDistancePopulation : public CalculateEuclideanDistancePopulation
{
public:
  std::vector<int> selected_residues_for_rmsd;
  std::vector<double> inter_distance_individual_1, inter_distance_individual_2;
  std::vector<std::vector<double> > inter_distances_per_ind;
  std::map<std::pair<int, int>, double > inter_dist_norm_max;
  CalculateEuclideanMarioPartialDistancePopulation() : CalculateEuclideanDistancePopulation() {
  }

  CalculateEuclideanMarioPartialDistancePopulation(FitFunctionPtr fitfunc);

  CalculateEuclideanMarioPartialDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius) ;

  void
  apply_to_population(const std::vector<Individual>& popul, std::vector<double>& shared_fitness, std::vector<double>& rmsd_to_native,std::vector<double>& distances_of_population, std::vector<std::vector<NeighStruct> >& neigh_per_ind  );

  double
  current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2);

  double
  distance_of_individual_partial_mario_euclidean(core::pose::PoseOP pose_ind ,const std::vector<Individual>& popul, double& shared_acc, std::vector<NeighStruct>& neighs, std::vector<double>& ind_distances);

  void
  build_inter_distances_of_straight();

  void
  build_inter_distances_of_population( const std::vector<Individual>& popul);
};




typedef boost::shared_ptr<CalculateDistancePopulation> CalculateDistancePopulationPtr;

#endif /* CALCULATERMSDDISTANCEPOPULATION_H */
