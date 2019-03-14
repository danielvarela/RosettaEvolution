#ifndef MOVER_DIFFERENTIALEVOLUTION_HPP
#define MOVER_DIFFERENTIALEVOLUTION_HPP

#include <boost/shared_ptr.hpp>

#include <boost/property_tree/ptree.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <string>

#include "DE_types.hpp"
#include "FitFunction.hpp"
#include "InitPopulation.hpp"
#include "PrintBestIndividual.hpp"
#include "CalculateRmsdDistancePopulation.hpp"
#include "LocalSearchIndividualMover.hpp"
#if(MPI_ENABLED)
#include "../../mpi_files/MasterRosettaCalculator.hpp"
#endif


class PopulComparer {
public:
  PopulComparer() {}

  bool operator() (Individual i, Individual j) {
    return ( i.score < j.score );
  }
};


class ConfigurationDE
{
public:
  int NP;
  double CR, F;
  double fit_rad;
  std::string prot;

  ConfigurationDE(int n, double c, double f, double rad, std::string p) : NP(n), CR(c), F(f), fit_rad(rad), prot(p) {
  }
};

class MoverDE
{
public:
  PrintBestIndividualPtr print_best_ind;
  CalculateDistancePopulationPtr calculate_distances_popul;
  boost::shared_ptr<PartialRMSDcalculator> calculator_partial_rmsd;
  bool use_print_class;
  bool new_best_found;
  Individual current_best_ind;
  int trial_sucess_n;
  double previous_best;
  PopulComparer popul_comparer;
  boost::property_tree::ptree app_options;
  std::vector<int> desired_gens; 
  int gen_count;
  int clean_inds_gen_limit;
  int local_search_popul_module;
  std::string id;
  std::string mutation_strategy;
  std::string crossover_strategy;
  std::string select_parents_strategy;
  std::vector<core::pose::PoseOP> popul_pdb;


#if(MPI_ENABLED)
  boost::shared_ptr<MasterRosettaCalculator> mpi_calculator;
#endif

  MoverDE();

  MoverDE(ConfigurationDE pt, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population);
  MoverDE( boost::property_tree::ptree options , FitFunctionPtr scfxn_in, std::vector<Individual> initial_population);

  void print_timestamp();

  void set_scfxn(FitFunctionPtr scfxn_in);

  void init_popul(std::vector<Individual>& popul);
  int get_best_rand_ind(int target, std::string select_parents_strategy);

  int get_single_point_crossover(std::string ss);

  virtual bool select_population(const std::vector<Individual>& trial_popul);

  virtual void print_generation_information(int gen_count, bool new_best_found);

  virtual void apply();

  void apply(std::vector<Individual>& popul_in, int NP_in, int Gmax_in);

  void reset_stat();

  bool update_stat(int i, Individual ind, std::vector<Individual> popul);

  virtual double score(Individual& ind);

  virtual  Individual sample_new_individual(int i, const Parents& p, std::vector<Individual> popul);

  int size();

  virtual void select_parents(int i, Parents& parent);


  virtual void print_fitness_population(int gen_count);
  virtual  void print_distances_population();
  void print_individual_information(int gen_count, bool new_best_found);

  Individual best_ind() {
    std::cout << "best " << best_idx << std::endl;
    return scfxn->convert(popul[best_idx]);
  }

  std::vector<Individual> popul;
  int NP;
  int Gmax;
  int D;
  double CR, F;

  double fit_rad;
  double avg_acc, best;
  int best_idx;
  int last_gen_best;
  double diff_avg_bef, diff_avg_aft;
  std::string prot;

  FitFunctionPtr scfxn;
  InitPopulationPtr init_popul_;
};

class SharedFitnessIndividual : public Individual {
public:
  boost::shared_ptr<Individual> ind;

  int index;
  double fit_shared;
  double rmsd;
  std::vector<NeighStruct> neighs;

  SharedFitnessIndividual(boost::shared_ptr<Individual> ind_input) :
    ind(ind_input) {
  }
};


class RescaledFitnessComparer {
public:
  RescaledFitnessComparer() {}

  bool operator() (SharedFitnessIndividual i, SharedFitnessIndividual j) {
    return ( i.fit_shared < j.fit_shared );
  }
};

class RescaledSharedFitnessComparer {
public:
  RescaledSharedFitnessComparer() {}

  bool operator() (SharedFitnessIndividual i, SharedFitnessIndividual j) {
    return ( i.fit_shared < j.fit_shared );
  }
};

class SharedMoverDE : public MoverDE
{
public:
  std::vector<double> shared_fitness, rmsd_to_native,   ind_inter_distance;
  std::vector<Individual> copy_popul_fitness_for_print;
  std::vector<double> copy_popul_shared_fitness_for_print;
  std::vector<std::vector<NeighStruct> > neighs_per_ind;
  RescaledSharedFitnessComparer cmp_shared_fitness;

  SharedMoverDE() : MoverDE() {}

  SharedMoverDE(ConfigurationDE pt, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population) : MoverDE(pt, scfxn_in, initial_population) {}
  SharedMoverDE(boost::property_tree::ptree options, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population) : MoverDE(options, scfxn_in, initial_population) {}

  bool select_population(const std::vector<Individual>& trial_popul);


  void print_fitness_population(int gen_count);

  void print_distances_population();
};


class HybridMoverDE : public MoverDE
{
public:
  std::vector<double> shared_fitness, rmsd_to_native,   ind_inter_distance;

  std::vector<std::vector<NeighStruct> > neighs_per_ind;
  RescaledFitnessComparer cmp_score;
  RescaledSharedFitnessComparer cmp_shared_fitness;
  boost::shared_ptr<LocalSearchIndividualMover> local_search;

  HybridMoverDE(ConfigurationDE pt, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in) : MoverDE(pt, scfxn_in, initial_population) {
    local_search = local_search_in;
  }
  HybridMoverDE(boost::property_tree::ptree options, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in) : MoverDE(options, scfxn_in, initial_population) {
    local_search = local_search_in;
  }

  void apply_local_search_at_population();

  void apply();

  void print_generation_information(int gen_count, bool new_best_found);
};


class SharedHybridMoverDE : public SharedMoverDE
{
public:
  std::vector<double> shared_fitness, rmsd_to_native,   ind_inter_distance;
  RescaledFitnessComparer cmp_score;
  RescaledSharedFitnessComparer cmp_shared_fitness;
  boost::shared_ptr<LocalSearchIndividualMover> local_search;

  SharedHybridMoverDE(ConfigurationDE pt, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in) : SharedMoverDE(pt, scfxn_in, initial_population) {
    local_search = local_search_in;
  }
  SharedHybridMoverDE(boost::property_tree::ptree options, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in) : SharedMoverDE(options, scfxn_in, initial_population) {
    local_search = local_search_in;
  }

  void apply_local_search_at_population();

  void apply();
  double score(Individual& ind);

  void local_search_apply(Individual& ind) {
    local_search->apply(ind);
  }

};

class CrowdingMoverDE : public SharedMoverDE
{
public:
  CrowdingMoverDE(ConfigurationDE pt, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population) : SharedMoverDE(pt, scfxn_in, initial_population) {}
  CrowdingMoverDE(boost::property_tree::ptree options, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population) : SharedMoverDE(options, scfxn_in, initial_population) {}

  bool select_population(const std::vector<Individual>& trial_popul);

  void print_fitness_population(int gen_count);
};


class CrowdingHybridMoverDE : public SharedHybridMoverDE
{
public:
  std::vector<double> shared_fitness, rmsd_to_native,   ind_inter_distance;
  boost::shared_ptr<LocalSearchIndividualMover> local_search;
  double current_cde_radius;

  CrowdingHybridMoverDE(ConfigurationDE pt, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in) : SharedHybridMoverDE(pt, scfxn_in, initial_population, local_search_in) {
    local_search = local_search_in;
  }
  CrowdingHybridMoverDE(boost::property_tree::ptree options, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in) : SharedHybridMoverDE(options, scfxn_in, initial_population, local_search_in) {
    local_search = local_search_in;
  }


  bool select_population(const std::vector<Individual>& trial_popul);
  void print_distances_population();

  void print_fitness_population(int gen_count);
};


class ResetOldIndsHybridDE : public HybridMoverDE
{
public:
  boost::shared_ptr<LocalSearchIndividualMover> local_search;
  boost::shared_ptr<LocalSearchIndividualMover> greedy_mover;
  std::vector<Individual> old_population;

  ResetOldIndsHybridDE(ConfigurationDE pt, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in, boost::shared_ptr<LocalSearchIndividualMover> greedy_mover_in ) : HybridMoverDE(pt, scfxn_in, initial_population, local_search_in) {
    local_search = local_search_in;
    greedy_mover = greedy_mover_in;
  }

  ResetOldIndsHybridDE(boost::property_tree::ptree options, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in, boost::shared_ptr<LocalSearchIndividualMover> greedy_mover_in) : HybridMoverDE(options, scfxn_in, initial_population, local_search_in) {
    greedy_mover = greedy_mover_in;
    local_search = local_search_in;
  }
  void apply();
  bool select_population(const std::vector<Individual>& trial_popul);
  void clean_old_inds();
};

class ResetOldIndsCrowdingHybridDE : public ResetOldIndsHybridDE
{
public:
  ResetOldIndsCrowdingHybridDE(ConfigurationDE pt, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in, boost::shared_ptr<LocalSearchIndividualMover> greedy_mover_in ) : ResetOldIndsHybridDE(pt, scfxn_in, initial_population, local_search_in, greedy_mover_in) {
  }
  ResetOldIndsCrowdingHybridDE(boost::property_tree::ptree options, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in, boost::shared_ptr<LocalSearchIndividualMover> greedy_mover_in) : ResetOldIndsHybridDE(options, scfxn_in, initial_population, local_search_in, greedy_mover_in) {
  }
  bool select_population(const std::vector<Individual>& trial_popul);
};


#if(MPI_ENABLED)

class MPICrowdingMoverDE : public CrowdingHybridMoverDE
{
public:

  boost::shared_ptr<LocalSearchIndividualMover> local_search;
  MPICrowdingMoverDE(ConfigurationDE pt, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in) : CrowdingHybridMoverDE(pt, scfxn_in, initial_population, local_search_in) {
    local_search = local_search_in;
    mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterRosettaCalculator(scfxn_in));
  }
  MPICrowdingMoverDE(boost::property_tree::ptree options, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in) : CrowdingHybridMoverDE(options, scfxn_in, initial_population, local_search_in) {
    local_search = local_search_in;
    mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterRosettaCalculator(scfxn_in));
  }

  double score(Individual& ind);
  void apply_local_search_at_population();
  void analysis_of_population(std::vector<Individual> popul, std::vector<core::pose::PoseOP> popul_pdb);

  void apply();
};

class MPISeedsMoverDE : public MPICrowdingMoverDE
{
public:

  boost::shared_ptr<LocalSearchIndividualMover> local_search;
  MPISeedsMoverDE(ConfigurationDE pt, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in) : MPICrowdingMoverDE(pt, scfxn_in, initial_population, local_search_in) {
    local_search = local_search_in;
    mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterRosettaCalculator(scfxn_in));
  }
  MPISeedsMoverDE(boost::property_tree::ptree options, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in) : MPICrowdingMoverDE(options, scfxn_in, initial_population, local_search_in) {
    local_search = local_search_in;
    mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterRosettaCalculator(scfxn_in));
  }

  std::vector<std::vector<Individual> > create_seeds(std::vector<Individual> popul, std::vector<core::pose::PoseOP> popul_pdb);
  void apply();
};


class MPIResetOldCrowdingHybridDE : public ResetOldIndsCrowdingHybridDE
{
public:
  boost::shared_ptr<LocalSearchIndividualMover> local_search;

  MPIResetOldCrowdingHybridDE(ConfigurationDE pt, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in, boost::shared_ptr<LocalSearchIndividualMover> greedy_search_in ) : ResetOldIndsCrowdingHybridDE(pt, scfxn_in, initial_population, local_search_in, greedy_search_in) {
    local_search = local_search_in;
    mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterRosettaCalculator(scfxn_in));
  }
  MPIResetOldCrowdingHybridDE(boost::property_tree::ptree options, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in, boost::shared_ptr<LocalSearchIndividualMover> greedy_search_in) : ResetOldIndsCrowdingHybridDE(options, scfxn_in, initial_population, local_search_in, greedy_search_in) {
    local_search = local_search_in;
    mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterRosettaCalculator(scfxn_in));
  }

  void apply();
  double score(Individual& ind);

  boost::shared_ptr<MasterRosettaCalculator> mpi_calculator;
};
#endif

// class ExtraSharedHybridMoverDE : public SharedHybridMoverDE
// {
// public:
//   std::vector<Individual> extra_popul;

//   ExtraSharedHybridMoverDE(FitFunctionPtr scfxn_in, InitPopulationPtr init_popul_in,  boost::shared_ptr<LocalSearchIndividualMover> local_search_in) : SharedHybridMoverDE(scfxn_in, init_popul_in, local_search_in) {
//   }

//   void apply();
// };


#endif // MOVER_DIFFERENTIALEVOLUTION_HPP
