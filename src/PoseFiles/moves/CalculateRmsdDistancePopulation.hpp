
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


class CalculateRmsdDistancePopulation
{
public:
  core::pose::PoseOP pose_ind, pose_other;
  FitFunctionPtr scorefxn;
  std::string ss;
  core::pose::PoseOP native_;
  boost::shared_ptr<PoseFunction> pfunc;
  double fit_radius = 3.0;

  CalculateRmsdDistancePopulation() {
  }

  CalculateRmsdDistancePopulation(const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius) ;

  class DistancesResult
  {
  public:
    DistancesResult() {}

    std::vector<double> shared_fitness, rmsd_to_native, distances_of_population;
    std::vector<int> neigh_per_ind;
  };

  std::vector<double>
  run_rmsd_to_native(const std::vector<Individual>& popul) ;

  DistancesResult
  run(const std::vector<Individual>& popul) ;

  void
  calculate_rmsd(const std::vector<Individual>& popul,    std::vector<double>& rmsd_to_native ) ;

  std::vector<int>
  find_nearest_parent(int ind_idx, Individual target, const std::vector<Individual>& popul) ;

  int
  find_nearest(Individual target, const std::vector<Individual>& popul) ;

  virtual void
  apply_to_population(const std::vector<Individual>& popul, std::vector<double>& shared_fitness, std::vector<double>& rmsd_to_native,std::vector<double>& distances_of_population, std::vector<int>& neigh_per_ind  ) ;

  virtual double rmsd_between_inds(Individual ind_left, Individual ind_right) ;

  double distance_of_individual(core::pose::PoseOP pose_ind,const std::vector<Individual>& popul, double& shared_acc, int& neighs, std::vector<double>& ind_distances) ;

  void print_distance_population( const std::vector<Individual>& popul ) ;
};

class CalculateEuclideanDistancePopulation : public CalculateRmsdDistancePopulation
{
public:
  CalculateEuclideanDistancePopulation() : CalculateRmsdDistancePopulation() {
  }

  CalculateEuclideanDistancePopulation(FitFunctionPtr fitfunc);

  CalculateEuclideanDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius) ;

  virtual void
  apply_to_population(const std::vector<Individual>& popul, std::vector<double>& shared_fitness, std::vector<double>& rmsd_to_native,std::vector<double>& distances_of_population, std::vector<int>& neigh_per_ind  ) ;

  virtual  double euclidean_two_individuals(const Individual& ind_left, const Individual& ind_right) ;

  double euclidean_distance(const std::vector<double>& v1, const std::vector<double>& v2) ;

  double rmsd_between_inds(Individual ind_left, Individual ind_right) ;

  double distance_of_individual_euclidean(Individual pose_ind,const std::vector<Individual>& popul, double& shared_acc, int& neighs, std::vector<double>& ind_distances) ;
};

class CalculateEuclideanLoopDistancePopulation : public CalculateEuclideanDistancePopulation
{
public:
  CalculateEuclideanLoopDistancePopulation() : CalculateEuclideanDistancePopulation() {
  }

  CalculateEuclideanLoopDistancePopulation(FitFunctionPtr fitfunc);

  CalculateEuclideanLoopDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius) ;

  double euclidean_two_individuals(const Individual& ind_left, const Individual& ind_right) ;

};


class CalculateEuclideanDiffAbsDistancePopulation : public CalculateEuclideanDistancePopulation
{
public:
  CalculateEuclideanDiffAbsDistancePopulation() : CalculateEuclideanDistancePopulation() {
  }

  CalculateEuclideanDiffAbsDistancePopulation(FitFunctionPtr fitfunc);

  CalculateEuclideanDiffAbsDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius) ;

  double euclidean_two_individuals(const Individual& ind_left, const Individual& ind_right) ;

};



class CalculateEuclideanMarioDistancePopulation : public CalculateEuclideanDistancePopulation
{
public:
  std::vector<int> selected_residues_for_rmsd;

  CalculateEuclideanMarioDistancePopulation() : CalculateEuclideanDistancePopulation() {
  }

  CalculateEuclideanMarioDistancePopulation(FitFunctionPtr fitfunc);

  CalculateEuclideanMarioDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius) ;

  double euclidean_two_individuals(const Individual& ind_left, const Individual& ind_right) ;

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
  apply_to_population(const std::vector<Individual>& popul, std::vector<double>& shared_fitness, std::vector<double>& rmsd_to_native,std::vector<double>& distances_of_population, std::vector<int>& neigh_per_ind  );

  double
distance_of_individual_partial_mario_euclidean(Individual pose_ind,const std::vector<Individual>& popul, double& shared_acc, int& neighs, std::vector<double>& ind_distances);

  double euclidean_two_individuals(const Individual& ind_left, const Individual& ind_right) ;
 void
  build_inter_distances_of_straight();

 void
  build_inter_distances_of_population( const std::vector<Individual>& popul);

};


class CalculateEuclideanRmsdWithoutSuperpDistancePopulation : public CalculateEuclideanDistancePopulation
{
public:
  CalculateEuclideanRmsdWithoutSuperpDistancePopulation() : CalculateEuclideanDistancePopulation() {
  }

  CalculateEuclideanRmsdWithoutSuperpDistancePopulation(FitFunctionPtr fitfunc);

  CalculateEuclideanRmsdWithoutSuperpDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius) ;

  double euclidean_two_individuals(const Individual& ind_left, const Individual& ind_right) ;

};





typedef boost::shared_ptr<CalculateRmsdDistancePopulation> CalculateRmsdDistancePopulationPtr;

#endif /* CALCULATERMSDDISTANCEPOPULATION_H */
