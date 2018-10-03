
#include "CalculateRmsdDistancePopulation.hpp"
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <vector>
#include <string>
#include <map>
#include <utility>

#define NORM_MULTIPLIED_FACTOR 0.5f

InterDistancesGraph::InterDistancesGraph() {
  for (double i = 0; i < 5; i = i + 0.1) {
    bin_values.push_back(i);
    occurrence_map[bin_values.back()] = 0;
  }
}

void
InterDistancesGraph::add_occurrence(double value) {
  bool inserted = false;
  for (int i = 1; i < bin_values.size() && !inserted; i++) {
    if (value < bin_values[i]) {
      occurrence_map[bin_values[i - 1]]++;
      inserted = true;
    }
  }
  if (!inserted) {
    occurrence_map[bin_values.back()]++;
  }
}

std::string
InterDistancesGraph::print_graph() {
  std::string result = "";
  for (std::map<double, int>::iterator i = occurrence_map.begin(); i != occurrence_map.end(); ++i) {
    result += std::to_string(i->first) + " " + std::to_string(i->second) + " ";
  }
  return result;
}

CalculateDistancePopulation::CalculateDistancePopulation(const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius)  {
  pose_ind = p->clone();
  pose_other = p->clone();
  scorefxn = sfxn;
  ss = ss_in;
  native_ = p->clone();
  fit_radius = radius;
  pfunc = boost::dynamic_pointer_cast<PoseFunction >(scorefxn);
}

std::vector<double>
CalculateDistancePopulation::run_rmsd_to_native(const std::vector<Individual>& popul) {
  std::vector<double> rmsd_to_native;

  build_pdb_population(popul, popul_pdb);

  for (int i = 0; i < popul_pdb.size(); i++) {
    rmsd_to_native.push_back( core::scoring::CA_rmsd(*popul_pdb[i], *native_));
  }

  return rmsd_to_native;
}

void
CalculateDistancePopulation::build_pdb_population(const std::vector<Individual>& population, std::vector<core::pose::PoseOP>& popul_pdb) {
    popul_pdb.resize(0);
    core::pose::PoseOP local_pose;
    for (int i = 0; i < population.size(); ++i) {
      local_pose = pose_ind->clone();
      pfunc->fill_pose(local_pose, population[i], ss);
      popul_pdb.push_back(local_pose);
    }
  }


CalculateDistancePopulation::DistancesResult
CalculateDistancePopulation::run(const std::vector<Individual>& popul) {
  CalculateDistancePopulation::DistancesResult result;
  std::vector<double> shared_fitness, rmsd_to_native, distances_of_population;
  std::vector<std::vector<NeighStruct> > neigh_per_ind;


  build_pdb_population(popul, popul_pdb);
  
  apply_to_population(popul, shared_fitness, rmsd_to_native, distances_of_population, neigh_per_ind);
  result.shared_fitness = shared_fitness;
  result.rmsd_to_native = rmsd_to_native;
  result.distances_of_population = distances_of_population;
  result.neigh_per_ind = neigh_per_ind;

  return result;
}

std::vector<int>
CalculateDistancePopulation::find_nearest_parent(int ind_idx, Individual target, const std::vector<Individual>& popul) {
  int popul_size = popul.size();
  core::pose::PoseOP pose_target = pose_ind->clone();
  pfunc->fill_pose(pose_target, target, ss);
  
  std::vector<core::pose::PoseOP> popul_pdb;
  build_pdb_population(popul, popul_pdb);
 
  int nearest_idx = 0;
  double best_found = 100000;
  std::map<int, double> individual_distances;
  for (int i = 0; i < popul_size; i++) {
    if (i != ind_idx) {
      double curr_score = current_distance_calculation(*popul_pdb[i], *pose_target);
      individual_distances[i] = curr_score;
      if (curr_score < best_found) {
	best_found = curr_score;
	nearest_idx = i;
      }
    }
  }
  // Declaring the type of Predicate that accepts 2 pairs and return a bool
  typedef std::function<bool(std::pair<int, double>, std::pair<int, double>)> Comparator;

  // Defining a lambda function to compare two pairs. It will compare two pairs using second fields
  Comparator compFunctor = [](std::pair<int, double> elem1 ,std::pair<int, double> elem2)
    {
      return elem1.second < elem2.second;
    };
  // Declaring a set that will store the pairs using above comparision logic
  std::set<std::pair<int, double>, Comparator> setOfInds(
							 individual_distances.begin(), individual_distances.end(), compFunctor);
  int top_individuals = 5;
  std::vector<int> solution;
  int i = 0;
  for (std::pair<int, double> element : setOfInds) {
    solution.push_back(element.second);
    i++;
    if (i == top_individuals) {
      break;
    }
  }
  return solution;
}
double  CalculateDistancePopulation::distance_of_individual(core::pose::PoseOP pose_ind,const std::vector<core::pose::PoseOP>& popul_pdb, double& shared_acc, std::vector<NeighStruct>& neighs, std::vector<double>& ind_distances) {
  double total_dist = 0.0, ind_dist = 0.0, sh_aux = 0.0;
  int popul_pdb_size = popul_pdb.size();
  neighs.resize(0);

  for (int j = 0; j < popul_pdb_size; j++) {
    ind_dist = current_distance_calculation(*popul_pdb[j], *pose_ind);
    ind_distances.push_back(ind_dist);
    if (ind_dist < fit_radius) {
      sh_aux += (1 - pow( (ind_dist / fit_radius) , 2 ));
      //neigh
      NeighStruct n(j, ind_dist);
      neighs.push_back(n);
    }
    total_dist += ind_dist;
  }

  if (total_dist > 0) {
    total_dist = total_dist / popul_pdb_size;
  } else {
    total_dist = 0.0;
  }

  shared_acc = sh_aux;

  return total_dist;
}



int
CalculateDistancePopulation::find_nearest(Individual target, const std::vector<Individual>& popul) {
  build_pdb_population(popul, popul_pdb);
 
  int popul_size = popul.size();
  core::pose::PoseOP pose_target = pose_ind->clone();
  pfunc->fill_pose(pose_target, target, ss);

  int nearest_idx = 0;
  double best_found = 100000;
  for (int i = 0; i < popul_size; i++) {
    double curr_dist = current_distance_calculation(*popul_pdb[i], *pose_target);
    if (curr_dist < best_found) {
      best_found = curr_dist;
      nearest_idx = i;
    }
  }
  return nearest_idx;
}


void
CalculateDistancePopulation::apply_to_population(const std::vector<Individual>& popul, std::vector<double>& shared_fitness, std::vector<double>& rmsd_to_native,std::vector<double>& distances_of_population, std::vector<std::vector<NeighStruct> >& neigh_per_ind  ) {
  double avg_distance_of_ind = 0.0;
  std::vector<double> ind_distance_avg;
  std::vector<double> shared_fitness_ind;
  std::vector<double> rmsd_to_native_vector;
  double dist_avg;
  int popul_size = popul_pdb.size();
  double sh;
  std::vector<double> ind_distances;
  for (int i = 0; i < popul_size; i++) {
    std::vector<NeighStruct> neighs;
    pose_ind = popul_pdb[i];
    rmsd_to_native_vector.push_back( core::scoring::CA_rmsd(*pose_ind, *native_));
    dist_avg = distance_of_individual(pose_ind, popul_pdb, sh, neighs, ind_distances);
    ind_distance_avg.push_back(dist_avg);
    shared_fitness_ind.push_back(sh);
    neigh_per_ind.push_back(neighs);
  }

  shared_fitness = shared_fitness_ind;
  rmsd_to_native = rmsd_to_native_vector;
  distances_of_population = ind_distances;
}

void
CalculateDistancePopulation::calculate_rmsd(const std::vector<Individual>& popul,    std::vector<double>& rmsd_to_native ) {
  std::vector<double> rmsd_to_native_vector;
  int popul_size = popul.size();

  for (int i = 0; i < popul_size; i++) {
    pfunc->fill_pose(pose_ind, popul[i], ss);
    rmsd_to_native_vector.push_back( core::scoring::CA_rmsd(*pose_ind, *native_));
  }

  rmsd_to_native = rmsd_to_native_vector;
}


double CalculateDistancePopulation::distance_between_inds(Individual ind_left, Individual ind_right) {
  core::pose::PoseOP pose_1 = pose_ind->clone(), pose_2 = pose_ind->clone();
  pfunc->fill_pose(pose_1, ind_left, ss);
  pfunc->fill_pose(pose_1, ind_right, ss);
  return current_distance_calculation(*pose_1, *pose_2);
}

double
CalculateRmsdDistancePopulation::current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2) {
  return core::scoring::CA_rmsd(pose_1, pose_2);
}

CalculateNativeDiffDistancePopulation::CalculateNativeDiffDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius)  : CalculateDistancePopulation(p, sfxn, ss_in, pscore, radius) {
}


double
CalculateNativeDiffDistancePopulation::current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2) {
  double right_native =  core::scoring::CA_rmsd(*native_, pose_1);
  double left_native =  core::scoring::CA_rmsd(*native_, pose_2);
  //  std::cout << "native_diff_abs " << std::abs(right_native - left_native) << std::endl;
  return std::abs(right_native - left_native);
}


void
CalculateDistancePopulation::print_distance_population( const std::vector<Individual>& popul ) {
  // std::vector<double> ind_distance_avg = apply_to_population(popul);
  // std::cout << "[DIST_POP] ";
  // for (int i = 0; i < ind_distance_avg.size(); i++) {
  //   std::cout << ind_distance_avg[i] << " ";
  // }
  // std::cout << std::endl;
}


CalculateEuclideanDistancePopulation::CalculateEuclideanDistancePopulation(FitFunctionPtr fitfunc) : CalculateDistancePopulation() {
  scorefxn = fitfunc;
}

CalculateEuclideanDistancePopulation::CalculateEuclideanDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius)  : CalculateDistancePopulation(p, sfxn, ss_in, pscore, radius) {
}

double CalculateEuclideanDistancePopulation::current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2) {
  double sum = 0.0;
  std::vector<double> vars_1, vars_2;
  for (int res = 1; res <= pose_1.total_residue(); ++res) {
    vars_1.push_back(pose_1.phi(res));
    vars_1.push_back(pose_1.psi(res));
    vars_2.push_back(pose_2.phi(res));
    vars_2.push_back(pose_2.psi(res));
  }
  
  double max_euc = 0;
  double sum_aux = 0;
  for (int i = 0; i < vars_1.size(); i++) {
    double diff = vars_1[i] - vars_2[i];
    if (diff > 180.0) {
      diff -= 360.0;
    }
    if (diff < -180.0) {
      diff += 360.0;
    }

    if (diff > max_euc) {
      max_euc = diff;
    }
    sum_aux += (std::abs(diff) / 180.0);
    sum += std::pow(diff / 180.0, 2);
  }
  double result = std::sqrt(sum);
  result = result / vars_1.size();
  result = result * 100;
  return result;
}


CalculateEuclideanLoopDistancePopulation::CalculateEuclideanLoopDistancePopulation(FitFunctionPtr fitfunc) : CalculateEuclideanDistancePopulation(fitfunc) {
}

CalculateEuclideanLoopDistancePopulation::CalculateEuclideanLoopDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius)  : CalculateEuclideanDistancePopulation(p, sfxn, ss_in, pscore, radius) {
}

double
CalculateEuclideanLoopDistancePopulation::current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2) {
  double sum = 0.0;
  std::vector<double> vars_1, vars_2;
  for (int res = 1; res <= pose_1.total_residue(); ++res) {
    vars_1.push_back(pose_1.phi(res));
    vars_1.push_back(pose_1.psi(res));
    vars_2.push_back(pose_2.phi(res));
    vars_2.push_back(pose_2.psi(res));
  }
  double max_euc = 0;
  double sum_aux = 0;
  int res_counter = 0;
  double total_euc = 0, total_euc_aux = 0;
  int total_res_counter = 0;
  for (int i = 0; i < vars_1.size(); i++) {
    double diff = vars_1[i] - vars_2[i];
    if (diff > 180.0) {
      diff -= 360.0;
    }
    if (diff < -180.0) {
      diff += 360.0;
    }

    total_euc_aux +=  (std::abs(diff) / 180.0);
    total_euc += std::pow(diff / 180.0, 2);
    if (ss[res_counter] == 'L') {
      sum_aux += (std::abs(diff) / 180.0);
      sum += std::pow(diff / 180.0, 2);
      total_res_counter++;
    }

    if (i % 2 != 0) {
      res_counter++;
    }
  }

  double result = std::sqrt(sum);

  //result = result / total_res_counter;
  result = result / vars_1.size();

  int multiplied_factor = 100;
  result = result * multiplied_factor;
  std::cout << "total_euc_aux " << (total_euc_aux / vars_1.size()) * multiplied_factor <<  " total " << (( std::sqrt(total_euc) / vars_1.size() ) * multiplied_factor) << " sum_aux " << (sum_aux / vars_1.size()) * multiplied_factor   << " only_loops " << result << std::endl;
  return result;
}

CalculateEuclideanDiffAbsDistancePopulation::CalculateEuclideanDiffAbsDistancePopulation(FitFunctionPtr fitfunc) : CalculateEuclideanDistancePopulation(fitfunc) {
}

CalculateEuclideanDiffAbsDistancePopulation::CalculateEuclideanDiffAbsDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius)  : CalculateEuclideanDistancePopulation(p, sfxn, ss_in, pscore, radius) {
}

double
CalculateEuclideanDiffAbsDistancePopulation::current_distance_calculation(core::pose::Pose &pose_1, core::pose::Pose &pose_2) {
  std::vector<double> vars_1, vars_2;
  for (int res = 1; res <= pose_1.total_residue(); ++res) {
    vars_1.push_back(pose_1.phi(res));
    vars_1.push_back(pose_1.psi(res));
    vars_2.push_back(pose_2.phi(res));
    vars_2.push_back(pose_2.psi(res));
  }

  return euclidean_two_vectors(vars_1, vars_2);
}

double
CalculateEuclideanDistancePopulation::euclidean_two_vectors(const std::vector<double> & vars_1, const std::vector<double>& vars_2) {
  double sum = 0.0;
  double max_euc = 0;
  double sum_aux = 0;
  int res_counter = 0;
  double total_euc = 0, total_euc_aux = 0;
  int total_res_counter = 0;
  int count_angles_diff = 0;
  int angles_at_loops = 0;
  for (int i = 0; i < vars_1.size(); i++) {
    double diff = vars_1[i] - vars_2[i];
    if (diff > 180.0) {
      diff -= 360.0;
    }
    if (diff < -180.0) {
      diff += 360.0;
    }

   total_euc_aux +=  (std::abs(diff) / 180.0);
    total_euc += std::pow(diff / 180.0, 2);

    if (ss[res_counter] == 'L') {
      angles_at_loops++;
      if (std::abs(diff) > 10.0) {
	count_angles_diff++;
      }

      sum_aux += (std::abs(diff) / 180.0);
      sum += std::pow(diff / 180.0, 2);
      total_res_counter++;
    }

    if (i % 2 != 0) {
      res_counter++;
    }
  }

  // core::pose::PoseOP pose_ind_1, pose_ind_2;
  // pose_ind_1 = native_->clone();
  // pose_ind_2 = native_->clone();
  // pfunc->fill_pose(pose_ind_1, ind_1, ss);
  // pfunc->fill_pose(pose_ind_2, ind_2, ss);
  // double rmsd = core::scoring::CA_rmsd(*pose_ind_1, *pose_ind_2);
  int multiplied_factor = 10;
  // double  result =  (sum_aux / vars_1.size() ) * multiplied_factor ;
  double  result =  (sum_aux / angles_at_loops  ) * multiplied_factor ;
  // std::cout << "[IMP] count_angles_diff " << count_angles_diff << " abs_diff " << result << " rmsd " << rmsd << std::endl;
  // double result = std::sqrt(sum);
  // //result = result / total_res_counter;
  // result = result / vars_1.size();
  // result = result * multiplied_factor;
  //  std::cout << "total_euc_aux " << (total_euc_aux / vars_1.size()) * multiplied_factor <<  " total " << (( std::sqrt(total_euc) / vars_1.size() ) * multiplied_factor) << " sum_aux " << (sum_aux / vars_1.size()) * multiplied_factor   << " only_loops " << result << std::endl;
  return result;
}

CalculateEuclideanMarioDistancePopulation::CalculateEuclideanMarioDistancePopulation(FitFunctionPtr fitfunc) : CalculateEuclideanDistancePopulation(fitfunc) {
}

CalculateEuclideanMarioDistancePopulation::CalculateEuclideanMarioDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius)  : CalculateEuclideanDistancePopulation(p, sfxn, ss_in, pscore, radius) {
  for (int i = 1; i < ss.size(); i++) {
    if ((ss[i] != 'L') && (ss[i -1] != ss[i])) {
      //std::cout << "ss at " << i << " val " << ss[i] << std::endl;
      selected_residues_for_rmsd.push_back(i + 1);
    }
    if ((ss[i] == 'L') && (ss[i -1] != ss[i])) {
      //      std::cout << "ss at " << i - 1 << " val " << ss[i - 1] << std::endl;
      selected_residues_for_rmsd.push_back(i);
    }
  }
}

double
CalculateEuclideanMarioDistancePopulation::current_distance_calculation(core::pose::Pose &pose_1, core::pose::Pose &pose_2) {
  double sum = 0.0;
  for ( int i = 0; i < selected_residues_for_rmsd.size(); ++i)  {
    core::Size res_num = core::Size(selected_residues_for_rmsd[i]);
    core::PointPosition calpha1_pos  = pose_1.residue(res_num).xyz("CA");
    core::PointPosition calpha2_pos = pose_2.residue(res_num).xyz("CA");
    sum += calpha1_pos.distance_squared(calpha2_pos);
 }
  sum = sum / selected_residues_for_rmsd.size();
  double result = std::sqrt(sum);
  return result;
}

CalculateEuclideanMarioPartialDistancePopulation::CalculateEuclideanMarioPartialDistancePopulation(FitFunctionPtr fitfunc) : CalculateEuclideanDistancePopulation(fitfunc) {
}

CalculateEuclideanMarioPartialDistancePopulation::CalculateEuclideanMarioPartialDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius)  : CalculateEuclideanDistancePopulation(p, sfxn, ss_in, pscore, radius) {
  bool init_found = false;
  bool end_found = false;
  int init_residue = 0, end_residue = 0;
  for (int i = 1; i < ss.size(); i++) {
    if ((ss[i] != 'L') && (ss[i -1] != ss[i])) {
      std::cout << "ss at " << i << " val " << ss[i] << std::endl;
      init_residue = i + 1;
      init_found = true;
    }
    if ((ss[i] == 'L') && (ss[i -1] != ss[i])) {
      std::cout << "ss at " << i - 1 << " val " << ss[i - 1] << std::endl;
      end_residue = i - 1;
      end_found = true;
    }
    if (init_found && end_found) {
      init_found = false;
      end_found = false;
      int half_residue = std::floor( (end_residue + init_residue) / 2 );
      selected_residues_for_rmsd.push_back( half_residue );
    }
  }

  // for (int i = 0; i < selected_residues_for_rmsd.size(); i++) {
  //   std::cout << "residue " << selected_residues_for_rmsd[i] << " " << ss[selected_residues_for_rmsd[i] - 1 ] << std::endl;
  // }

  // int paa;
  // std::cin >> paa;

}

void
CalculateEuclideanMarioPartialDistancePopulation::build_inter_distances_of_straight() {
  core::pose::PoseOP straight = pose_ind->clone();

  for (int i = 1; i <= straight->total_residue(); i++) {
    straight->set_phi(i, 180.0);
    straight->set_psi(i, 180.0);
  }
  int half_res = straight->total_residue() / 2;
  core::PointPosition calpha1_pos_aux  = straight->residue(1).xyz("CA");
  core::PointPosition calpha2_pos_aux = straight->residue(half_res).xyz("CA");
  double max_dist = calpha1_pos_aux.distance(calpha2_pos_aux);

    for ( int i = 0; i < selected_residues_for_rmsd.size(); ++i)  {
      for (int j  = i + 1; j < selected_residues_for_rmsd.size(); j++) {
	core::Size res_num_1 = core::Size(selected_residues_for_rmsd[i]);
	core::Size res_num_2 = core::Size(selected_residues_for_rmsd[j]);
	core::PointPosition calpha1_pos  = straight->residue(res_num_1).xyz("CA");
	core::PointPosition calpha2_pos = straight->residue(res_num_2).xyz("CA");
	double dist = calpha1_pos.distance(calpha2_pos);
	// falta normalizar esta distancia dividiendo por al distancia maxima ( la misma dist para la proteina estirada)
	//inter_distance_individual_1.push_back( std::sqrt( std::pow(dist, 2) / 2)  );

	inter_dist_norm_max[std::make_pair<int, int>(res_num_1, res_num_2)] = max_dist ;
      }
    }
}

void
CalculateEuclideanMarioPartialDistancePopulation::apply_to_population(const std::vector<Individual>& popul, std::vector<double>& shared_fitness, std::vector<double>& rmsd_to_native,std::vector<double>& distances_of_population, std::vector<std::vector<NeighStruct> >& neigh_per_ind  ) {
  double avg_distance_of_ind = 0.0;
  std::vector<double> ind_distance_avg;
  std::vector<double> shared_fitness_ind;
  std::vector<double> rmsd_to_native_vector;
  double dist_avg;
  int popul_size = popul.size();
  double sh;
  std::vector<double> ind_distances;

  // build inter_distances_per_ind;
  build_inter_distances_of_straight();
  build_inter_distances_of_population(popul);

  for (int i = 0; i < popul_size; i++) {
    std::vector<NeighStruct> neighs;
    rmsd_to_native_vector.push_back( core::scoring::CA_rmsd(*popul_pdb[i], *native_));

    inter_distance_individual_1 = inter_distances_per_ind[i];
    dist_avg = distance_of_individual_partial_mario_euclidean(popul_pdb[i], popul, sh, neighs, ind_distances );
    ind_distance_avg.push_back(dist_avg);
    shared_fitness_ind.push_back(sh);
    neigh_per_ind.push_back(neighs);
  }

  shared_fitness = shared_fitness_ind;
  rmsd_to_native = rmsd_to_native_vector;
  distances_of_population = ind_distances;
}

double
CalculateEuclideanMarioPartialDistancePopulation::distance_of_individual_partial_mario_euclidean(core::pose::PoseOP pose_ind,const std::vector<Individual>& popul, double& shared_acc,std::vector<NeighStruct>& neighs, std::vector<double>& ind_distances) {
  double euc = 0.0, ind_euc = 0.0, sh_aux = 0.0;
  //    double fit_radius = 0.001;
  //double fit_radius = 0.001;
  int popul_size = popul.size();
  neighs.resize(0);
  double min_inter_distance = 1000, max_inter_distance = -1000;
  for (int j = 0; j < popul_size; j++) {
    inter_distance_individual_2 = inter_distances_per_ind[j];
    ind_euc = current_distance_calculation(*pose_ind, *popul_pdb[j]);

    if (min_inter_distance > ind_euc) {
      min_inter_distance = ind_euc;
    }
    if (max_inter_distance < ind_euc) {
      max_inter_distance = ind_euc;
    }

    ind_distances.push_back(ind_euc);

    if (ind_euc < fit_radius) {
      double sh_acc =  (1 -  (ind_euc / fit_radius));
      sh_aux += sh_acc;
      NeighStruct n(j, ind_euc);
      neighs.push_back(n);
    }
    euc += ind_euc;
  }
  // std::cout << "!!!! min_inter_distance " << min_inter_distance << " max_inter_distance " << max_inter_distance << std::endl;

  if (euc > 0) {
    euc = euc / popul_size;
  } else {
    euc = 0.0;
  }

  shared_acc = sh_aux;

  return euc;
}

void
CalculateEuclideanMarioPartialDistancePopulation::build_inter_distances_of_population( const std::vector<Individual>& popul) {
  inter_distances_per_ind.resize(0);
  for (int i = 0; i < popul.size(); i++) {
    core::pose::PoseOP pose_ind_1, pose_ind_2;
    pose_ind_1 = pose_ind->clone();
    pfunc->fill_pose(pose_ind_1, popul[i], ss);
    std::vector<double> inter_distance_individual;
    build_inter_distance_for_an_individual(pose_ind_1, inter_distance_individual);
    inter_distances_per_ind.push_back(inter_distance_individual);
    // std::cout << "inter_ind " << i << " : ";
    // for (int d = 0; d < inter_distance_individual_1.size(); d++) {
    //   std::cout << inter_distance_individual_1[d] << " , ";
    // }
    // std::cout << std::endl;
  }
}

void
CalculateEuclideanMarioPartialDistancePopulation::build_inter_distance_for_an_individual(core::pose::PoseOP pose_ind, std::vector<double>& inter_distance_individual) {
  inter_distance_individual.resize(0);
  for ( int i = 0; i < selected_residues_for_rmsd.size(); ++i)  {
    for (int j  = i + 1; j < selected_residues_for_rmsd.size(); j++) {
      core::Size res_num_1 = core::Size(selected_residues_for_rmsd[i]);
      core::Size res_num_2 = core::Size(selected_residues_for_rmsd[j]);
      core::PointPosition calpha1_pos  = pose_ind->residue(res_num_1).xyz("CA");
      core::PointPosition calpha2_pos = pose_ind->residue(res_num_2).xyz("CA");
      double dist = calpha1_pos.distance(calpha2_pos);
      // normalizar esta distancia dividiendo por al distancia maxima ( la misma dist para la proteina estirada)
      //inter_distance_individual.push_back( std::sqrt( std::pow(dist, 2) / 2)  );
      if (inter_dist_norm_max[std::pair<int, int>(res_num_1, res_num_2)] == 0) {
	inter_distance_individual.push_back(0);
      } else {
	double normalized_dist = dist / ( NORM_MULTIPLIED_FACTOR * inter_dist_norm_max[std::pair<int, int>(res_num_1, res_num_2)] );
	//	  std::cout << " dist " << dist << " norm " << normalized_dist << " max " << inter_dist_norm_max[std::pair<int, int>(res_num_1, res_num_2)] << std::endl;
	inter_distance_individual.push_back( normalized_dist );
      }
    }
  }
}

int
CalculateEuclideanMarioPartialDistancePopulation::find_nearest(Individual target, const std::vector<Individual>& popul) {
  build_pdb_population(popul, popul_pdb);
  build_inter_distances_of_population(popul);

  int popul_size = popul.size();
  core::pose::PoseOP pose_target = pose_ind->clone();
  pfunc->fill_pose(pose_target, target, ss);

  std::vector<double> inter_distance_target;
  build_inter_distance_for_an_individual(pose_target, inter_distance_target);

  int nearest_idx = 0;
  double best_found = 100000;
  inter_distance_individual_1 = inter_distance_target;
  // std::cout << "find_nearest " << std::endl;
  for (int i = 0; i < popul_size; i++) {
    inter_distance_individual_2 = inter_distances_per_ind[i];
    double curr_dist = current_distance_calculation(*popul_pdb[i], *pose_target);
    // std::cout << i << " : " << curr_dist << std::endl;
    if (curr_dist < best_found) {
      best_found = curr_dist;
      nearest_idx = i;
    }
  }

  // std::cout<< "found at " << nearest_idx << " dist " << best_found << std::endl;

  return nearest_idx;
}

double
CalculateEuclideanMarioPartialDistancePopulation::current_distance_calculation(core::pose::Pose& pose_1, core::pose::Pose& pose_2) {
  double sum = 0.0;
  // std::cout << "number of distances " << inter_distance_individual_1.size() << std::endl;
  for ( int i = 0; i < inter_distance_individual_1.size(); ++i)  {
    sum += std::pow( inter_distance_individual_1[i] - inter_distance_individual_2[i], 2);
    // std::cout << inter_distance_individual_1[i] << " - " << inter_distance_individual_2[i] << std::endl;
    // std::cout << "acum: " << sum << std::endl;
  }

  double result = std::sqrt(sum);
  result = result / inter_distance_individual_1.size();
  // std::cout << "final result after sqrt and norm " << result << std::endl;
  // int p;
  // std::cin >> p;

  return result;
}


