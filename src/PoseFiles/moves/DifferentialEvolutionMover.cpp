#include "DifferentialEvolutionMover.hpp"

#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <string>
#include <omp.h>

#include <stdlib.h>
#include <boost/make_shared.hpp>

#define BOOST_DATE_TIME_NO_LIB

// Includes
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/local_time_adjustor.hpp>
#include <boost/date_time/c_local_time_adjustor.hpp>
#include <boost/optional/optional.hpp>

#include "../../helper_apps/FileToPDBPrinter.hpp"

MoverDE::MoverDE() {
  srand(time(NULL));
  use_print_class = false;
  //scfxn = boost::shared_ptr<FitFunction>(new DeJongFunction());
  //scfxn = boost::shared_ptr<FitFunction>(new HimmelblauFunction());
  scfxn = boost::shared_ptr<FitFunction>(new ToyFunction());
  init_popul_ = boost::shared_ptr<InitPopulation>(new InitPopulation(scfxn));

  NP = 100;
  Gmax = 100;
  D = scfxn->D();
  CR = 0.9;
  F = 0.5;

  avg_acc = 0;
  best = 0;
  best_idx = 0;

  // init_popul(popul);
}

MoverDE::MoverDE(ConfigurationDE pt, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population) {
  srand(time(NULL));
  use_print_class = false;
  scfxn = scfxn_in;
  //  init_popul_ = init_popul_in;

  NP = pt.NP;
  Gmax = 10000;
  D = scfxn->D();

  CR = pt.CR;
  F = pt.F;
  prot = pt.prot;

  avg_acc = 1;
  best = 100000;
  best_idx = 0;
  last_gen_best = 0;
  fit_rad = pt.fit_rad;
  //  init_popul(popul);
  popul = initial_population;
  current_best_ind = popul[0];
  std::cout << "[CONF] ";
  std::cout << " NP " << NP ;
  std::cout << " CR " << CR ;
  std::cout << " F " << F ;
  std::cout << " rad " << pt.fit_rad;
  print_timestamp();
  std::cout << std::endl;
}


MoverDE::MoverDE(boost::property_tree::ptree options, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population) {
  srand(time(NULL));
  app_options = options;
  
  use_print_class = false;
  scfxn = scfxn_in;
  //  init_popul_ = init_popul_in;

  NP = app_options.get<int>("DE.NP");
  Gmax = 10000;
  D = scfxn->D();

  CR = app_options.get<double>("DE.CR");
  F = app_options.get<double>("DE.F");
  prot = app_options.get<std::string>("Protocol.prot");

  boost::property_tree::ptree::const_assoc_iterator it = app_options.find("PrintPopul.gens");
  if(app_options.not_found() == it) {
    std::string desired_gens_opt = app_options.get<std::string>("PrintPopul.gens");
    boost::char_separator<char> sep(",");
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    tokenizer tokens(desired_gens_opt,sep);
    for (tokenizer::iterator it = tokens.begin(); it != tokens.end(); ++it) {
      desired_gens.push_back(std::stoi(*it));
    }
  }

  avg_acc = 1;
  best = 100000;
  best_idx = 0;
  last_gen_best = 0;
  fit_rad = app_options.get<double>("Extra.fitrad");
  //  init_popul(popul);
  popul = initial_population;
  current_best_ind = popul[0];
  std::cout << "[CONF] ";
  std::cout << " NP " << NP ;
  std::cout << " CR " << CR ;
  std::cout << " F " << F ;
  std::cout << " rad " << fit_rad;
  print_timestamp();
  std::cout << std::endl;
}

void MoverDE::print_timestamp() {
  // Code
  using boost::posix_time::ptime;
  using boost::posix_time::second_clock;
  using boost::posix_time::to_simple_string;
  using boost::gregorian::day_clock;

  ptime todayUtc(day_clock::universal_day(), second_clock::universal_time().time_of_day());
  std::cout << " time " << to_simple_string(todayUtc) << std::endl;

}

void MoverDE::set_scfxn(FitFunctionPtr scfxn_in) {
  scfxn = scfxn_in;
}

void MoverDE::init_popul(std::vector<Individual>& popul) {
  init_popul_->apply(popul, NP, size());
}

void MoverDE::apply() {
  gen_count = 0;
  Parents parents;
  new_best_found = false;
  std::cout << "init DE for function " << scfxn->name() << std::endl;
  std::vector<Individual> trial_popul;

  while (gen_count < Gmax) {
    reset_stat();
    trial_popul.resize(0);

     // #pragma omp parallel
    // #pragma omp for
    for (int i = 0; i < NP; ++i) {
      score(popul[i]);
      select_parents(i, parents);

      Individual ind = sample_new_individual(i, parents, popul);
      trial_popul.push_back(ind);
    }

    new_best_found = select_population(trial_popul);

    print_generation_information(gen_count, new_best_found);
    // std::ofstream myfile;
    // std::string filename = "./.matplotlib_files/output_popul_gen_"+std::to_string(gen_count)+".txt";
    // myfile.open (filename);
    // for (int i = 0; i < popul.size(); i++) {
    //   Individual ind_converted = scfxn->convert(popul[i]);
    //   myfile << std::setprecision(4) << ind_converted.vars[0] << " " << ind_converted.vars[1] << " " << ind_converted.score << std::endl;
    // }
    // myfile.close();
    // std::cout << "End of gen " << gen_count << std::endl;
    // int p;
    // std::cin >> p;
    ++gen_count;

  }
}

void HybridMoverDE::apply() {
  gen_count = 0;
  Parents parents;
  new_best_found = false;
  std::cout << "init Hybrid DE for function " << scfxn->name() << std::endl;
  std::vector<Individual> trial_popul;
  double global_min_rmsd, global_max_rmsd, global_avg_rmsd;
  while (gen_count < Gmax) {
    reset_stat();
    local_search->reset_stats();
    trial_popul.resize(0);
    global_min_rmsd = 100000;
    global_max_rmsd = -100000;
    global_avg_rmsd = 0;

    //#pragma omp parallel for private(local_search)
    for (int i = 0; i < NP; ++i) {
      //std::cout<<"threads="<<omp_get_num_threads()<<std::endl;

      // apparently it disturbs the right score of the population
      local_search->apply(popul[i]);
    }

    for (int i = 0; i < NP; ++i) {
      select_parents(i, parents);

      Individual ind = sample_new_individual(i, parents, popul);
      if (CR > 0.5) {
	double new_rmsd = calculate_distances_popul->distance_between_inds(ind, popul[parents.x1]);
	if (new_rmsd > global_max_rmsd) {
	  global_max_rmsd = new_rmsd;
	}
	if (new_rmsd < global_min_rmsd) {
	  global_min_rmsd = new_rmsd;
	}
	global_avg_rmsd += new_rmsd;
      } else {
	double new_rmsd = calculate_distances_popul->distance_between_inds(ind, popul[i]);
	std::cout << "rmsd parent " << i << " " << new_rmsd << std::endl;
      }

      trial_popul.push_back(ind);
    }

    std::cout << "[RMSD_TRIALS] GEN " << gen_count << " { " << global_min_rmsd << " , " << global_max_rmsd << " } " << global_avg_rmsd / NP << std::endl;

    std::cout << "[TRIAL_0] " << "total " << trial_popul[0].score << " " << print_best_ind->print_energies(trial_popul[0]) << std::endl;
    new_best_found = select_population(trial_popul);

    print_generation_information(gen_count, new_best_found);
    ++gen_count;
  }
}


bool MoverDE::select_population(const std::vector<Individual>& trial_popul) {
  bool new_best_found = false;
  std::vector<Individual> large_popul = trial_popul;
  large_popul.insert(large_popul.end(), popul.begin(), popul.end());

  std::sort(large_popul.begin(), large_popul.end(), popul_comparer);
  large_popul.resize(popul.size());

  if (large_popul[0].score < popul[0].score) {
    new_best_found = true;
  }

  popul = large_popul;

  // default version of selection
  //   trial_sucess_n = 0;
//   for (int i = 0; i < NP; i++) {
//     if (trial_popul[i].score < popul[i].score) {
//       popul[i] = trial_popul[i];
//       trial_sucess_n++;
//     }

//
//   }

//   return new_best_found;
  new_best_found = update_stat(0, trial_popul[0], popul) || new_best_found;
  return new_best_found;
}

bool CrowdingHybridMoverDE::select_population(const std::vector<Individual>& trial_popul) {
  trial_sucess_n = 0;
  // Crowding DE: Rene Thomsen Algorithm
  // 1. Each trial can be replace only its nearest individual if the energy is better

  bool new_best_found = false;
  double previous_best = popul[best_idx].score;
  for (int i = 0; i < trial_popul.size(); i++) {
    int nearest_ind = calculate_distances_popul->find_nearest(trial_popul[i], popul);
    if ( trial_popul[i].score < popul[nearest_ind].score ) {
      popul[nearest_ind] = trial_popul[i];
      new_best_found = true;
      trial_sucess_n++;
    }
  }

  avg_acc = 0;
  for (int i = 0; i < popul.size(); i++) {
    avg_acc += popul[i].score;
    if (popul[i].score < popul[best_idx].score) {
      best_idx = i;
      best = popul[best_idx].score;
      if (best < previous_best) {
	new_best_found = true;
      }
    }
  }

  return new_best_found;
}


void CrowdingMoverDE::print_fitness_population(int gen_count) {
  //if ((gen_count < 25) || (gen_count % 50 == 0)) {
  if (true) {
    std::cout << "[POP] ";
    for (int i = 0; i < popul.size(); i++) {
      std::cout << (-1 * SCORE_ERROR_FIXED) + popul[i].score << " ";
    }
    std::cout << std::endl;
  }
}


bool CrowdingMoverDE::select_population(const std::vector<Individual>& trial_popul) {
  trial_sucess_n = 0;
  // Crowding DE: Rene Thomsen Algorithm
  // 1. Each trial can be replace only its nearest individual if the energy is better

  bool new_best_found = false;
  double previous_best = popul[best_idx].score;
  for (int i = 0; i < trial_popul.size(); i++) {
    int nearest_ind = calculate_distances_popul->find_nearest(trial_popul[i], popul);
    double previous_score = popul[nearest_ind].score;
    if ( trial_popul[i].score < popul[nearest_ind].score ) {
      std::cout << "nearest ind " << nearest_ind << " ind " << i << "(" << calculate_distances_popul->distance_between_inds(trial_popul[i], popul[nearest_ind])   << ") : " << trial_popul[i].score << " " << previous_score << std::endl;
      popul[nearest_ind] = trial_popul[i];
      int p;
      std::cin >> p;
      new_best_found = true;
      trial_sucess_n++;
    }
  }

  avg_acc = 0;
  for (int i = 0; i < popul.size(); i++) {
    avg_acc += popul[i].score;
    if (popul[i].score < popul[best_idx].score) {
      best_idx = i;
      best = popul[best_idx].score;
      if (best < previous_best) {
	new_best_found = true;
      }
    }
  }

  return new_best_found;
}

bool SharedMoverDE::select_population(const std::vector<Individual>& trial_popul) {
  bool new_best_found = false;
  trial_sucess_n = 0;
  // Shared DE: Rene Thomsen Algorithm
  // 1. Join trial and current popul into a large population
  // 2. calculate fitness sharing for each individual of this population
  // 3. Sort the population according to its fitness sharing
  // 4. Select the best half of the population to continue in the evolution
  //     4.1. If the selected population doen't contains the best individual, then change it for the last worst individual according to its shared fitness
  std::vector<Individual> large_popul = trial_popul;
  if (desired_gens.begin() != desired_gens.end()) {
    if (std::find(desired_gens.begin(), desired_gens.end(), gen_count) !=
	desired_gens.end()) {
      FileToPDBPrinter::population_to_file(app_options, scfxn->name(),gen_count, trial_popul, FileToPDBPrinter::PrintPopulationOption::print_trial_population);
    }
  }

  large_popul.insert(large_popul.end(), popul.begin(), popul.end());

  double using_sf = false;

  using_sf = true;
  CalculateDistancePopulation::DistancesResult result =
    calculate_distances_popul->run(large_popul);

  shared_fitness = result.shared_fitness;
  rmsd_to_native = result.rmsd_to_native;
  neighs_per_ind = result.neigh_per_ind;
  ind_inter_distance = result.distances_of_population;

  std::vector<SharedFitnessIndividual> shared_popul;
  int large_best_idx = 0;
  double large_best_value = 100000000000000;
  for (int i = 0; i < large_popul.size(); i++) {
    SharedFitnessIndividual sf_ind = SharedFitnessIndividual(boost::make_shared<Individual>(large_popul[i]) );
    sf_ind.index = i;
    double current_score = large_popul[i].score;
    if (current_score > 0) {
      sf_ind.fit_shared = (current_score) *  shared_fitness[i] ;
    } else {
      sf_ind.fit_shared = (current_score) /  shared_fitness[i] ;
    }
    sf_ind.score = current_score;
    sf_ind.rmsd = rmsd_to_native[i];
    sf_ind.neighs = neighs_per_ind[i];

    shared_popul.push_back(sf_ind);
    if (current_score < large_best_value) {
      large_best_idx = i;
      large_best_value = current_score;
    }
  }

  new_best_found = (large_best_idx < NP);
  best = large_popul[large_best_idx].score;

  std::sort(shared_popul.begin(), shared_popul.end(), cmp_shared_fitness);

  bool best_already_included = false;
  avg_acc = 0;
  copy_popul_fitness_for_print.resize(large_popul.size());
  copy_popul_shared_fitness_for_print.resize(large_popul.size());
  for (int i = 0; i < shared_popul.size(); i++) {
    int include_index = shared_popul[i].index;
    if (i < NP) {
      if (include_index < NP) trial_sucess_n++;
      popul[i] = large_popul[include_index];
      //    popul[i].score = shared_popul[i].ind->score;
      shared_fitness[i] = shared_popul[i].fit_shared;
      avg_acc += popul[i].score;
      rmsd_to_native[i] = shared_popul[i].rmsd;
      if (include_index == large_best_idx) {
	best_already_included = true;
	best_idx = i;
      }
    }

    copy_popul_fitness_for_print[i] = large_popul[include_index];
    copy_popul_shared_fitness_for_print[i] = shared_popul[i].fit_shared;
    neighs_per_ind[i] = shared_popul[i].neighs;
  }

  if (!best_already_included) {
    trial_sucess_n++;
    popul[popul.size() - 1] = large_popul[large_best_idx];
    //    popul[popul.size() - 1].score = shared_popul[large_best_idx].ind->score;
    shared_fitness[popul.size() - 1] = shared_popul[large_best_idx].fit_shared;
    rmsd_to_native[popul.size() - 1] = shared_popul[large_best_idx].rmsd;
    best_idx = popul.size() - 1;
  }

  shared_fitness.resize(popul.size());
  rmsd_to_native.resize(popul.size());

  return new_best_found;
}

void MoverDE::print_generation_information(int gen_count, bool new_best_found) {

  std::cout << "[GEN] " << gen_count << " " << (-1 * SCORE_ERROR_FIXED) + best << " " << (-1 * SCORE_ERROR_FIXED) + (avg_acc / NP) << " at " << best_idx << " ";
  //std::cout << "[GEN] " << gen_count << " " << best << " " << avg_acc / NP << " at " << best_idx << " ";
  std::cout << std::endl;

  std::cout << scfxn->print_stats() << std::endl;

  if (scfxn->name() == "first De Jong") {
      print_fitness_population(gen_count);

      Individual ind_converted = scfxn->convert(popul[best_idx]);
      print_distances_population();
      std::cout << "[BEST_IND] " << popul[best_idx].score << " ";
      std::cout << "vars ";
      for (int i = 0; i < scfxn->D(); i++) {
	std::cout << 	ind_converted.vars[i] << " ";
      }
      std::cout << std::endl;
  }

  if ((getenv("PRINT_ALL") != NULL) && (scfxn->name() != "first De Jong")) {
    if (( std::string(getenv("PRINT_ALL")) == "true") && (gen_count % 1 == 0)) {
      print_fitness_population(gen_count);
      print_distances_population();
      print_individual_information(gen_count, new_best_found);

      std::cout << "[TRIAL_SUC] " << trial_sucess_n << std::endl;

      if (desired_gens.begin() != desired_gens.end()) {
          if (std::find(desired_gens.begin(), desired_gens.end(), gen_count) !=
              desired_gens.end()) {
            FileToPDBPrinter::population_to_file(app_options, scfxn->name(),gen_count, popul);
          }
        }


    }
  }
}

void HybridMoverDE::print_generation_information(int gen_count, bool new_best_found) {
  std::cout << "[GEN] " << gen_count << " " << (-1 * SCORE_ERROR_FIXED) + best << " " << (-1 * SCORE_ERROR_FIXED) + (avg_acc / NP) << " at " << best_idx << " ";
  std::cout << std::endl;

  std::cout << scfxn->print_stats() << std::endl;
  std::cout << local_search->print_stats() << std::endl;
  std::cout << "[TRIAL_SUC] " << trial_sucess_n << std::endl;

  if (getenv("PRINT_ALL") != NULL) {
    if (( std::string(getenv("PRINT_ALL")) == "true")) {
      print_fitness_population(gen_count);
      print_distances_population();
      print_individual_information(gen_count, new_best_found);
    }
  }
}



double MoverDE::score(Individual& ind) {
  return scfxn->score(ind);
}

Individual MoverDE::sample_new_individual(int i, const Parents& p, std::vector<Individual> popul) {
  Individual trial(size());

  int n_res = 0;
  for (int k = 0; k < D; ++k) {
    if (k > 1 && k % 2 != 0 ) {
      n_res++;
    }
    double jrand = URAND();
    //if ( ( jrand < CR ) && (popul[i].ss[n_res] != 'H') && (popul[i].ss[n_res] != 'E') ) {
    if ( ( jrand < CR ) ) {
      double diff =  popul[p.x2].vars[k] - popul[p.x3].vars[k];
      diff_avg_bef += std::abs(diff);
      while (diff > 1) { diff -= 1; }
      while (diff < -1) { diff += 1; }
      diff_avg_aft += std::abs(diff);
      trial.vars[k] = popul[p.x1].vars[k] + F* (diff);

      //           std::cout << "n_res " << n_res << " " << popul[i].ss[n_res] <<  std::endl;
      if (trial.vars[k] > DE_LIMIT) {
      	trial.vars[k] = DE_LIMIT - std::abs(trial.vars[k] - DE_LIMIT);
      }
      if (trial.vars[k] < (-1*DE_LIMIT)) {
      	trial.vars[k] = (-1) + std::abs(trial.vars[k] - (-1*DE_LIMIT));
      }
      
    } else {
      trial.vars[k] = popul[i].vars[k];
    }
  }

  trial.omega = popul[p.x1].omega;
  trial.ss = popul[p.x1].ss;

  score(trial);
  return trial;
}

int MoverDE::size() { return D; }

void MoverDE::select_parents(int i, Parents& parent) {
  do {
    parent.x1 = rand() % NP;
  } while (parent.x1 == i);


  //parent.x1 = best_idx;
  do {
    parent.x2 = rand() % NP;
  } while ((parent.x2 ==i) || (parent.x2 == parent.x1));

  do {
    parent.x3 = rand() % NP;
  } while ((parent.x3 == i) || (parent.x3 == parent.x2) || (parent.x3 == parent.x1));
}


void CrowdingMoverDE::select_parents(int i, Parents& parent) {
  do {
    parent.x1 = rand() % NP;
  } while (parent.x1 == i);
  do {
    parent.x2 = rand() % NP;
  } while ((parent.x2 ==i) || (parent.x2 == parent.x1));

  do {
    parent.x3 = rand() % NP;
  } while ((parent.x3 == i) || (parent.x3 == parent.x2) || (parent.x3 == parent.x1));
}


void CrowdingHybridMoverDE::select_parents(int i, Parents& parent) {

  // std::vector<int> nearest_inds = calculate_distances_popul->find_nearest_parent(i, popul[i], popul);

  do {
    parent.x1 = rand() % NP;
  } while (parent.x1 == i);

  // parent.x1 = best_idx;
  // parent.x1 = nearest_inds[0];
  // parent.x1 = nearest_inds[1];
  // parent.x1 = nearest_inds[2];

  do {
    parent.x2 = rand() % NP;
  } while ((parent.x2 ==i) || (parent.x2 == parent.x1));

  do {
    parent.x3 = rand() % NP;
  } while ((parent.x3 == i) || (parent.x3 == parent.x2) || (parent.x3 == parent.x1));
}



void MoverDE::print_fitness_population(int gen_count) {
  //if ((gen_count < 25) || (gen_count % 50 == 0)) {
  if (true) {
    std::cout << "[POP] ";
    for (int i = 0; i < NP; i++) {
      std::cout << (-1 * SCORE_ERROR_FIXED) + popul[i].score << " ";
    }
    std::cout << std::endl;
  }
}

void SharedMoverDE::print_fitness_population(int gen_count) {
  //if ((gen_count < 25) || (gen_count % 50 == 0)) {
  if (true) {
    std::cout << "[POP] ";
    for (int i = 0; i < copy_popul_fitness_for_print.size(); i++) {
      std::cout << (-1 * SCORE_ERROR_FIXED) + copy_popul_fitness_for_print[i].score << " ";
    }
    std::cout << std::endl;
  }
}


void SharedMoverDE::print_distances_population() {
  //std::vector<double> ind_distance_avg = shared_fitness;
  std::vector<double> ind_distance_avg = copy_popul_shared_fitness_for_print;


  CalculateDistancePopulation::DistancesResult result =
                                        calculate_distances_popul->run(popul);

  int large_NP = copy_popul_fitness_for_print.size();

  //solo ocurre con CrowdingHybridMoverDE
  if (ind_distance_avg.size() < NP) {
    ind_distance_avg = result.distances_of_population;
  }

  std::vector<double> rmsd_native = result.rmsd_to_native;
  std::vector<double> ind_inter_distance = result.distances_of_population;

  if (scfxn->name() != "first De Jong") {
    std::cout << "[RMSD_NATIVE] ";
    for (int i = 0; i < NP; i++) {
      std::cout << rmsd_native[i] << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "[DIST_POP] ";
  std::cout << ind_distance_avg[best_idx]  + (-1 * SCORE_ERROR_FIXED) << " ";
  for (int i = 0; i < large_NP; i++) {
    std::cout << ind_distance_avg[i] + (-1 * SCORE_ERROR_FIXED) << " ";
  }
  std::cout << std::endl;

  InterDistancesGraph inter_distances_calc;
  /* inter distances graphic */
  double min_inter_distance = 1000, max_inter_distance = -1000;
  for (int i = 0; i < ind_inter_distance.size(); i++) {
    // if (min_inter_distance > ind_inter_distance[i]) {
    //   min_inter_distance = ind_inter_distance[i];
    // }
    // if (max_inter_distance < ind_inter_distance[i]) {
    //   max_inter_distance =ind_inter_distance[i];
    // }
    //    std::cout << ind_inter_distance[i] << std::endl;
    inter_distances_calc.add_occurrence(ind_inter_distance[i]);
  }
  //  std::cout << "min_inter_distance " << min_inter_distance << " max_inter_distance " << max_inter_distance << std::endl;
  std::cout << "[INTER] ";
  std::cout << inter_distances_calc.print_graph() << std::endl;
  /* END inter distances graphic */

  if (neighs_per_ind.size() > 0) {
    std::cout << "[NEIGH] ";

    for (int i = 0; i < large_NP; i++) {
      std::cout << neighs_per_ind[i].size() << " ";
    }
    std::cout << std::endl;
  } else {
    std::cout << "[NEIGH] ";
    for (int i = 0; i < NP; i++) {
      std::cout << result.neigh_per_ind[i].size() << " ";
    }
    std::cout << std::endl;
  }
}


void MoverDE::print_distances_population() {
  // CalculateDistancePopulation::DistancesResult result =
  //                                       calculate_distances_popul->run(popul);

  // std::vector<double> ind_distance_avg = result.distances_of_population;
 // std::vector<double> rmsd_native = result.rmsd_to_native;
  std::vector<double> rmsd_native = calculate_distances_popul->run_rmsd_to_native(popul);
  // std::vector<double> ind_inter_distance = result.distances_of_population;

  std::cout << "[DIST_POP] ";
  std::cout << " ";
  //for (int i = 0; i < ind_distance_avg.size(); i++) {
  for (int i = 0; i < 50; i++) {
    std::cout << "1 ";
  }
  std::cout << std::endl;

  InterDistancesGraph inter_distances_calc;

  //for (int i = 0; i < ind_inter_distance.size(); i++) {
    for (int i = 0; i < 10; i++) {
            inter_distances_calc.add_occurrence(0);
    }

    std::cout << "[INTER] ";
    std::cout << inter_distances_calc.print_graph() << std::endl;


  std::cout << "[RMSD_NATIVE] ";
  for (int i = 0; i < NP; i++) {
    std::cout << rmsd_native[i] << " ";
  }
  std::cout << std::endl;
}

void MoverDE::print_individual_information(int gen_count, bool new_best_found) {
  int ind_index = best_idx;
  if (gen_count == 0) {
    previous_best = popul[best_idx].score;
  }

  if (previous_best > popul[best_idx].score) {
    new_best_found = true;
    previous_best = popul[best_idx].score;
  }

  if (use_print_class) {
    if ((new_best_found) || (gen_count == 0)) {
      double perc_of_diff = 0;
      for (int i = 0; i < D; i++) {
	if (std::abs(popul[ind_index].vars[i] - current_best_ind.vars[i]) > 0.013) {
	  perc_of_diff++;
	}
      }
      if (perc_of_diff > 0) {
	perc_of_diff = (perc_of_diff * 100)/D;
      }
      current_best_ind = popul[ind_index];

      std::string url_best = print_best_ind->print(popul[ind_index]);
      std::cout << "[IND] " << gen_count << " " << (-1 * SCORE_ERROR_FIXED) + popul[ind_index].score << " " << perc_of_diff << " " << url_best << std::endl;
      std::cout << "[ENERGY] " << print_best_ind->print_energies(popul[ind_index]) << std::endl;
    }
  }


}

void MoverDE::reset_stat() {
  avg_acc = 0;
  // best = 100000;
  //best_idx = 0;
  new_best_found = false;
  scfxn->reset_stats();
}

bool MoverDE::update_stat(int i, Individual ind, std::vector<Individual> popul) {
  avg_acc += popul[i].score;
  bool result = false;
  if (popul[i].score < best) {
    best_idx = i;
    best = popul[i].score;
    result = true;
  }
  return result;
}

void SharedHybridMoverDE::apply() {
  gen_count = 0;
  Parents parents;
  new_best_found = false;
  std::cout << "init Shared Hybrid DE for function " << scfxn->name() << std::endl;
  std::vector<Individual> trial_popul;
  double global_min_rmsd, global_max_rmsd, global_avg_rmsd;
  while (gen_count < Gmax) {
    reset_stat();
    local_search->reset_stats();
    trial_popul.resize(0);
    global_min_rmsd = 100000;
    global_max_rmsd = -100000;
    global_avg_rmsd = 0;

    double total_diff_avg_bef=0, total_diff_avg_aft = 0;
    diff_avg_bef=0; diff_avg_aft = 0;

    //    LocalSearchIndividualMover prt_local_search = *local_search.get();

    //    #pragma omp parallel for shared(prt_local_search)
    for (int i = 0; i < NP; ++i) {
      //      std::cout<<"threads="<<omp_get_num_threads()<<std::endl;
      // apparently it disturbs the right score of the population
      //prt_local_search.apply(popul[i]);
       local_search->apply(popul[i]);
    }

    for (int i = 0; i < NP; ++i) {
      select_parents(i, parents);

      Individual ind = sample_new_individual(i, parents, popul);
      diff_avg_bef = (diff_avg_bef / D);
      diff_avg_aft = (diff_avg_aft / D);
      total_diff_avg_bef += diff_avg_bef;
      total_diff_avg_aft += diff_avg_aft;
      if (CR > 0.5) {
	double new_rmsd = calculate_distances_popul->distance_between_inds(ind, popul[parents.x1]);
	//	std::cout << "rmsd parent " << i << " " << new_rmsd << std::endl;
	if (new_rmsd > global_max_rmsd) {
	  global_max_rmsd = new_rmsd;
	}
	if (new_rmsd < global_min_rmsd) {
	  global_min_rmsd = new_rmsd;
	}
	global_avg_rmsd += new_rmsd;
      } else {
	double new_rmsd = calculate_distances_popul->distance_between_inds(ind, popul[i]);
	std::cout << "rmsd parent " << i << " " << new_rmsd << std::endl;
      }

      trial_popul.push_back(ind);
    }

    std::cout << "[DIFF_AVG] " << (total_diff_avg_bef / NP) << " "  << (total_diff_avg_aft / NP) << std::endl;

    std::cout << "[RMSD_TRIALS] GEN " << gen_count << " { " << global_min_rmsd << " , " << global_max_rmsd << " } " << global_avg_rmsd / NP << std::endl;

    std::cout << "[TRIAL_0] " << "total " << trial_popul[0].score << " " << print_best_ind->print_energies(trial_popul[0]) << std::endl;
    new_best_found = select_population(trial_popul);

    print_generation_information(gen_count, new_best_found);
    ++gen_count;
  }
}


// otra forma para sample individual
/*      if (trial.vars[k] > DE_LIMIT) {
	trial.vars[k] = 1.0 - (URAND() * ( 1.0 * 0.05));
      }
      if (trial.vars[k] < (-1*DE_LIMIT)) {
	trial.vars[k] = -1.0 + (URAND() * ( -1.0 * 0.05));
      }


 */
