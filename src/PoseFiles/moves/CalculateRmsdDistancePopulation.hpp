
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

  enum distances_enum {
    dist_error, partial_rmsd, rmsd, rmsd_native_diff, euclidean, euclidean_loop, euclidean_diff_abs, euclidean_partial_mario, euclidean_mario, mario_first_last
  };
  std::map<std::string, distances_enum> distances_map;


  boost::shared_ptr<CalculateDistancePopulation>
  use_distances_strategy(std::string option);

  CalculateDistancePopulation() {
    init_distances_map();
  }

  void init_distances_map() {
    distances_map["error"] = dist_error;
    distances_map["rmsd"] = rmsd;
    distances_map["partial_rmsd"] = partial_rmsd;
    distances_map["rmsd_native_diff"] = rmsd_native_diff;
    distances_map["euclidean"] = euclidean;
    distances_map["euclidean_loop"] = euclidean_loop;
    distances_map["euclidean_diff_abs"] = euclidean_diff_abs;
    distances_map["euclidean_mario"] = euclidean_mario;
    distances_map["euclidean_partial_mario"] = euclidean_partial_mario;
    distances_map["mario_first_last"] = mario_first_last;
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

  virtual int
  find_nearest(Individual target, const std::vector<Individual>& popul, int target_index = -1, std::vector<core::pose::PoseOP> popul_pdb = std::vector<core::pose::PoseOP>());

  virtual std::vector<int>
  find_top_nearest(Individual target, const std::vector<Individual>& popul, int target_index, std::vector<core::pose::PoseOP> popul_pdb);



  std::vector<int>
  find_nearest_parent(int ind_idx, Individual target, const std::vector<Individual>& popul);

  double distance_between_inds(Individual ind_left, Individual ind_right);

  virtual double
  current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2) {
    return 10000343003;
  }

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
  build_inter_distance_for_an_individual(core::pose::PoseOP pose_ind_1, std::vector<double>& inter_distance_individual_1);

  int
  find_nearest(Individual target, const std::vector<Individual>& popul, int target_index = -1, std::vector<core::pose::PoseOP> popul_pdb = std::vector<core::pose::PoseOP>() );

  std::vector<int>
  find_top_nearest(Individual target, const std::vector<Individual>& popul, int target_index, std::vector<core::pose::PoseOP> popul_pdb);
 double
  individual_distance_value(Individual target, Individual reference);

  double
  pose_distance_calculation(core::pose::PoseOP pose1, core::pose::PoseOP pose2);


  void
  build_inter_distances_of_population( const std::vector<Individual>& popul);
};

class CalculateEuclideanMarioFirstLast : public CalculateEuclideanMarioPartialDistancePopulation  {
  CalculateEuclideanMarioFirstLast() : CalculateEuclideanMarioPartialDistancePopulation() {
  }

  CalculateEuclideanMarioFirstLast(FitFunctionPtr fitfunc) : CalculateEuclideanMarioPartialDistancePopulation(fitfunc)  {}

  CalculateEuclideanMarioFirstLast (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius) :  CalculateEuclideanMarioPartialDistancePopulation(p, sfxn, ss_in, pscore, radius)   {}

  double current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2);
};

class PartialRMSDcalculator: public CalculateDistancePopulation
{
public:
  PartialRMSDcalculator() : CalculateDistancePopulation() {
  }

  PartialRMSDcalculator(const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius) : CalculateDistancePopulation(p, sfxn, ss_in, pscore, radius) {}

  std::vector<std::pair<int, int> > selected_regions;

  double current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2){
    if (partial_distance(pose_1, pose_2).size() > 0) {
      return partial_distance(pose_1, pose_2)[0];
    } else {
      return 0;
    }

  }

  std::vector<double>
  partial_distance(core::pose::Pose &pose_1, core::pose::Pose &pose_2);

  std::vector<std::vector<double> >
  partial_rmsd_of_population(const std::vector<Individual>& popul);

  std::vector<std::pair<int, int> >
  make_selected_regions();

  std::string
  make_head(const std::vector<Individual>& popul);

};

typedef boost::shared_ptr<CalculateDistancePopulation> CalculateDistancePopulationPtr;

class PrintRMSDvsMarioAnalisys
{
public:
  class PrintRow
  {
  public:
    int i;
    int j;
    double rmsd;
    double euc;

    PrintRow() : i(0), j(0), rmsd(0), euc(0) {}
    PrintRow(int i_, int j_, double rmsd_, double euc_ ) : i(i_), j(j_), rmsd(rmsd_), euc(euc_) {
    }
    std::string print() {
      return std::string("[" + std::to_string(rmsd) + "," + std::to_string(euc) + "]" );
    }
  };

  CalculateDistancePopulationPtr calculate_partial_rmsd, calculate_partial_euc;
  std::vector<PrintRow> printable;

  PrintRMSDvsMarioAnalisys( CalculateDistancePopulationPtr calculate_partial_rmsd_, CalculateDistancePopulationPtr calculate_partial_euc_  ) : calculate_partial_rmsd(calculate_partial_rmsd_), calculate_partial_euc(calculate_partial_euc_) {
  }

  void apply(std::vector<Individual> popul) {
    std::vector<core::pose::PoseOP> popul_pdb;
    popul_pdb.resize(popul.size());
    boost::shared_ptr<CalculateEuclideanMarioPartialDistancePopulation> ffxn_partial_euc;
    ffxn_partial_euc = boost::shared_ptr<CalculateEuclideanMarioPartialDistancePopulation>( new CalculateEuclideanMarioPartialDistancePopulation( calculate_partial_rmsd->pose_ind, calculate_partial_rmsd->scorefxn, calculate_partial_rmsd->ss, NULL, calculate_partial_rmsd->fit_radius));
    ffxn_partial_euc->build_pdb_population(popul,popul_pdb);
    calculate_partial_rmsd->build_pdb_population(popul, popul_pdb);
    ffxn_partial_euc->build_inter_distances_of_population(popul);
    double euc = 0, rmsd = 0;
    printable.resize(0);
    for (int i = 0; i < popul.size() / 2 - 1; i++) {
      for (int j = i + 1; j < popul.size() / 2; j++) {
	euc = ffxn_partial_euc->individual_distance_value(popul[i], popul[j]);
	rmsd = calculate_partial_rmsd->current_distance_calculation(*popul_pdb[i], *popul_pdb[j] );
	PrintRow r(i, j, rmsd, euc);
	printable.push_back(r);
      }
    }

    std::cout << "[MARIO_PAIRS] ";
    for (int i = 0; i < printable.size(); i++) {
      std::cout << std::string(printable[i].print() + " ; ");
    }
    std::cout << std::endl;
  }
};




#endif /* CALCULATERMSDDISTANCEPOPULATION_H */
