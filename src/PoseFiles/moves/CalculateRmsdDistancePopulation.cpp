
#include "CalculateRmsdDistancePopulation.hpp"
#include <vector>
#include <string>
#include <map>

InterDistancesGraph::InterDistancesGraph() {
  for (double i = 0; i < 15; i = i + 0.33) {
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

CalculateRmsdDistancePopulation::CalculateRmsdDistancePopulation(const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius)  {
  pose_ind = p->clone();
  pose_other = p->clone();
  scorefxn = sfxn;
  ss = ss_in;
  native_ = p->clone();
  fit_radius = radius;
  pfunc = boost::dynamic_pointer_cast<PoseFunction >(scorefxn);
}


std::vector<double>
CalculateRmsdDistancePopulation::run_rmsd_to_native(const std::vector<Individual>& popul) {
  std::vector<double> rmsd_to_native;

  for (int i = 0; i < popul.size(); i++) {
    pfunc->fill_pose(pose_ind, popul[i], ss);
    rmsd_to_native.push_back( core::scoring::CA_rmsd(*pose_ind, *native_));
  }

  return rmsd_to_native;
}

CalculateRmsdDistancePopulation::DistancesResult
CalculateRmsdDistancePopulation::run(const std::vector<Individual>& popul) {
  CalculateRmsdDistancePopulation::DistancesResult result;
  std::vector<double> shared_fitness, rmsd_to_native, distances_of_population;
  std::vector<int> neigh_per_ind;

  apply_to_population(popul, shared_fitness, rmsd_to_native, distances_of_population, neigh_per_ind);
  result.shared_fitness = shared_fitness;
  result.rmsd_to_native = rmsd_to_native;
  result.distances_of_population = distances_of_population;
  result.neigh_per_ind = neigh_per_ind;

  return result;
}

void
CalculateRmsdDistancePopulation::calculate_rmsd(const std::vector<Individual>& popul,    std::vector<double>& rmsd_to_native ) {
  std::vector<double> rmsd_to_native_vector;
  int popul_size = popul.size();

  for (int i = 0; i < popul_size; i++) {
    pfunc->fill_pose(pose_ind, popul[i], ss);
    rmsd_to_native_vector.push_back( core::scoring::CA_rmsd(*pose_ind, *native_));
  }

  rmsd_to_native = rmsd_to_native_vector;
}


std::vector<int>
CalculateRmsdDistancePopulation::find_nearest_parent(int ind_idx, Individual target, const std::vector<Individual>& popul) {
  int popul_size = popul.size();
  core::pose::PoseOP pose_target = pose_ind->clone();
  pfunc->fill_pose(pose_target, target, ss);

  int nearest_idx = 0;
  double best_found = 100000;
  std::map<int, double> individual_distances;
  for (int i = 0; i < popul_size; i++) {
    if (i != ind_idx) {
      pfunc->fill_pose(pose_ind, popul[i], ss);
      double curr_score = core::scoring::CA_rmsd(*pose_ind, *pose_target) ;
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

int
CalculateRmsdDistancePopulation::find_nearest(Individual target, const std::vector<Individual>& popul) {
  int popul_size = popul.size();
  core::pose::PoseOP pose_target = pose_ind->clone();
  pfunc->fill_pose(pose_target, target, ss);

  int nearest_idx = 0;
  double best_found = 100000;
  for (int i = 0; i < popul_size; i++) {
    pfunc->fill_pose(pose_ind, popul[i], ss);
    double curr_score = core::scoring::CA_rmsd(*pose_ind, *pose_target) ;
    if (curr_score < best_found) {
      best_found = curr_score;
      nearest_idx = i;
    }
  }

  return nearest_idx;
}


void
CalculateRmsdDistancePopulation::apply_to_population(const std::vector<Individual>& popul, std::vector<double>& shared_fitness, std::vector<double>& rmsd_to_native,std::vector<double>& distances_of_population, std::vector<int>& neigh_per_ind  ) {
  double avg_distance_of_ind = 0.0;
  std::vector<double> ind_distance_avg;
  std::vector<double> shared_fitness_ind;
  std::vector<double> rmsd_to_native_vector;
  double rmsd_avg;
  int popul_size = popul.size();
  double sh;
  int neighs;


  std::vector<double> ind_distances;
  for (int i = 0; i < popul_size; i++) {
    pfunc->fill_pose(pose_ind, popul[i], ss);
    rmsd_to_native_vector.push_back( core::scoring::CA_rmsd(*pose_ind, *native_));
    rmsd_avg = distance_of_individual(pose_ind, popul, sh, neighs, ind_distances);
    ind_distance_avg.push_back(rmsd_avg);
    shared_fitness_ind.push_back(sh);
    neigh_per_ind.push_back(neighs);
  }

  shared_fitness = shared_fitness_ind;
  rmsd_to_native = rmsd_to_native_vector;
  distances_of_population = ind_distances;
}

double  CalculateRmsdDistancePopulation::rmsd_between_inds(Individual ind_left, Individual ind_right) {
  core::pose::PoseOP pose_ind_right = pose_ind->clone();
  pfunc->fill_pose(pose_ind, ind_left, ss);
  pfunc->fill_pose(pose_ind_right, ind_right, ss);
  return core::scoring::CA_rmsd(*pose_ind, *pose_ind_right);
}

double  CalculateRmsdDistancePopulation::distance_of_individual(core::pose::PoseOP pose_ind,const std::vector<Individual>& popul, double& shared_acc, int& neighs, std::vector<double>& ind_distances) {
  double rmsd = 0.0, ind_rmsd = 0.0, sh_aux = 0.0;



  int popul_size = popul.size();
  int neigh = 0;

  for (int j = 0; j < popul_size; j++) {
    pfunc->fill_pose(pose_other, popul[j], ss);
    ind_rmsd = core::scoring::CA_rmsd(*pose_other, *pose_ind);
    ind_distances.push_back(ind_rmsd);
    if (ind_rmsd < fit_radius) {
      sh_aux += (1 - pow( (ind_rmsd / fit_radius) , 2 ));
      neigh++;
    }
    rmsd += ind_rmsd;
  }

  if (rmsd > 0) {
    rmsd = rmsd / popul_size;
  } else {
    rmsd = 0.0;
  }


  shared_acc = sh_aux;
  neighs = neigh;

  return rmsd;
}

void
CalculateRmsdDistancePopulation::print_distance_population( const std::vector<Individual>& popul ) {
  // std::vector<double> ind_distance_avg = apply_to_population(popul);
  // std::cout << "[DIST_POP] ";
  // for (int i = 0; i < ind_distance_avg.size(); i++) {
  //   std::cout << ind_distance_avg[i] << " ";
  // }
  // std::cout << std::endl;
}


CalculateEuclideanDistancePopulation::CalculateEuclideanDistancePopulation(FitFunctionPtr fitfunc) : CalculateRmsdDistancePopulation() {
  scorefxn = fitfunc;
}

CalculateEuclideanDistancePopulation::CalculateEuclideanDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius)  : CalculateRmsdDistancePopulation(p, sfxn, ss_in, pscore, radius) {
}

void
CalculateEuclideanDistancePopulation::apply_to_population(const std::vector<Individual>& popul, std::vector<double>& shared_fitness, std::vector<double>& rmsd_to_native,std::vector<double>& distances_of_population, std::vector<int>& neigh_per_ind  ) {
  double avg_distance_of_ind = 0.0;
  std::vector<double> ind_distance_avg;
  std::vector<double> shared_fitness_ind;
  std::vector<double> rmsd_to_native_vector;
  double rmsd_avg;
  int popul_size = popul.size();
  double sh;
  int neighs;

  std::vector<double> ind_distances;
  for (int i = 0; i < popul_size; i++) {
    pfunc->fill_pose(pose_ind, popul[i], ss);
    rmsd_to_native_vector.push_back( core::scoring::CA_rmsd(*pose_ind, *native_));
    rmsd_avg = distance_of_individual_euclidean(popul[i], popul, sh, neighs, ind_distances );
    ind_distance_avg.push_back(rmsd_avg);
    shared_fitness_ind.push_back(sh);
    neigh_per_ind.push_back(neighs);
  }

  shared_fitness = shared_fitness_ind;
  rmsd_to_native = rmsd_to_native_vector;
  distances_of_population = ind_distances;
}

double CalculateEuclideanDistancePopulation::euclidean_two_individuals(const Individual& ind_left, const Individual& ind_right) {
  // COMENTAR ESTAS LINEAS
  double sum = 0.0;
  Individual ind_1 = scorefxn->convert(ind_left);
  Individual ind_2 = scorefxn->convert(ind_right);


  double max_euc = 0;
  double sum_aux = 0;
  for (int i = 0; i < ind_1.vars.size(); i++) {
    double diff = ind_1.vars[i] - ind_2.vars[i];
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

  result = result / ind_1.vars.size();

  result = result * 100;

  //    std::cout << "result acc " << sum_aux / ind_1.vars.size() << " euclidean " << result << " max_diff " << max_euc / 180.0 << std::endl;

  //    result = result * 10;

  return result;

  // Individual ind_1 = scorefxn->convert(ind_left);
  // Individual ind_2 = scorefxn->convert(ind_right);
  // return euclidean_distance(ind_1.vars, ind_2.vars);

  // ACORDARSE DE MODIFICAR ESTA LINEA
  //return euclidean_distance(ind_left.vars, ind_right.vars);
}

double CalculateEuclideanDistancePopulation::euclidean_distance(const std::vector<double>& v1, const std::vector<double>& v2) {
  // if vector size is the same
  if (v1.size() == v2.size())
    {
      double sumVal = 0.0;

      // sum the squared difference
      for (int i = 0;i < v1.size();i++)
	{
	  double diff = (double) (v1[i] - v2[i]);
	  sumVal += (diff*diff);
	}

      // if sum is bigger than zero
      if (sumVal > 0)
	{
	  // return the square root of the sum as the Euclidean distance
	  return sqrt(sumVal);
	}
      else // all differences were zero, so report 0.0 as Euclidean distance        {
	return 0.0;
    }
  else
    {
      throw std::invalid_argument( "Vector lengths differ in Euclidean distance calculation");
    }

}

double CalculateEuclideanDistancePopulation::rmsd_between_inds(Individual ind_left, Individual ind_right) {
  return euclidean_two_individuals(ind_left, ind_right);
}

double CalculateEuclideanDistancePopulation::distance_of_individual_euclidean(Individual pose_ind,const std::vector<Individual>& popul, double& shared_acc, int& neighs, std::vector<double>& ind_distances) {
  double euc = 0.0, ind_euc = 0.0, sh_aux = 0.0;
  //    double fit_radius = 0.001;
  //double fit_radius = 0.001;
  int popul_size = popul.size();
  int neigh = 0;

  for (int j = 0; j < popul_size; j++) {
    ind_euc = euclidean_two_individuals(pose_ind, popul[j]);
    ind_distances.push_back(ind_euc);
    if (ind_euc < fit_radius) {
      double sh_acc =  (1 -  (ind_euc / fit_radius));
      sh_aux += sh_acc;
      neigh++;
    }
    euc += ind_euc;
  }

  if (euc > 0) {
    euc = euc / popul_size;
  } else {
    euc = 0.0;
  }

  shared_acc = sh_aux;
  neighs = neigh;

  return euc;
}


CalculateEuclideanLoopDistancePopulation::CalculateEuclideanLoopDistancePopulation(FitFunctionPtr fitfunc) : CalculateEuclideanDistancePopulation() {
  scorefxn = fitfunc;
}

CalculateEuclideanLoopDistancePopulation::CalculateEuclideanLoopDistancePopulation (const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore, double radius)  : CalculateEuclideanDistancePopulation(p, sfxn, ss_in, pscore, radius) {
}

double
CalculateEuclideanLoopDistancePopulation::euclidean_two_individuals(const Individual& ind_left, const Individual& ind_right) {
  // COMENTAR ESTAS LINEAS
  double sum = 0.0;
  Individual ind_1 = scorefxn->convert(ind_left);
  Individual ind_2 = scorefxn->convert(ind_right);

  double max_euc = 0;
  double sum_aux = 0;
  int res_counter = 0;
  for (int i = 0; i < ind_1.vars.size(); i++) {
    if (ss[res_counter] == 'L') {
      double diff = ind_1.vars[i] - ind_2.vars[i];
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
    if (i % 2 == 0) {
      res_counter++;
    }
  }

  double result = std::sqrt(sum);

  result = result / ind_1.vars.size();

  result = result * 100;
  return result;
}
