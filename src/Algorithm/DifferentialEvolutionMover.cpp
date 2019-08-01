#include "DifferentialEvolutionMover.hpp"



#include <vector>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <string>
#include <omp.h>
#include <chrono>
#include <utility>
#include <map>
#include <boost/pointer_cast.hpp>
#include <boost/shared_ptr.hpp>

#include <stdlib.h>
#include <boost/make_shared.hpp>

#define BOOST_DATE_TIME_NO_LIB

// Includes
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/local_time_adjustor.hpp>
#include <boost/date_time/c_local_time_adjustor.hpp>
#include <boost/optional/optional.hpp>

#include "../Extra/FileToPDBPrinter.hpp"

MoverDE::MoverDE() {
  srand(time(NULL));
  use_print_class = false;
  //scfxn = boost::shared_ptr<FitFunction>(new DeJongFunction());
  //scfxn = boost::shared_ptr<FitFunction>(new HimmelblauFunction());
  scfxn = boost::shared_ptr<FitFunction>(new ToyFunction());
  init_popul_ = boost::shared_ptr<InitPopulation>(new InitPopulation(scfxn));

  id = std::to_string( static_cast<int>(std::abs(rand())) );

  NP = 100;
  Gmax = 100;
  D = scfxn->D();
  CR = 0.9;
  F = 0.5;

  avg_acc = 0;
  best = 0;
  best_idx = 0;
  mutation_strategy = std::string("default");
  clean_inds_gen_limit = 1000000;
  local_search_popul_module = 10;
  // init_popul(popul);
}

MoverDE::MoverDE(ConfigurationDE pt, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population) {
  srand(time(NULL));
  use_print_class = false;
  scfxn = scfxn_in;
  id = std::to_string( static_cast<int>(std::abs(rand())) );

#if(MPI_ENABLED)
  mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterRosettaCalculator(scfxn_in));
  std::string score_name = scfxn->name();
  std::string score_func = std::string(score_name.substr(score_name.rfind(" ") + 1));
  mpi_calculator->stage_name = score_func;
#endif

  mutation_strategy = std::string("default");
  //  init_popul_ = init_popul_in;

  NP = pt.NP;
  Gmax = 10000;
  D = scfxn->D();

  CR = pt.CR;
  F = pt.F;
  prot = pt.prot;
  local_search_popul_module = 10;

  avg_acc = 1;
  best = 100000;
  best_idx = 0;
  last_gen_best = 0;
  fit_rad = pt.fit_rad;
  //  init_popul(popul);
  popul = initial_population;
  current_best_ind = popul[0];
  best = popul[0].score;
  std::cout << "[CONF] ";
  std::cout << " NP " << NP ;
  std::cout << " CR " << CR ;
  std::cout << " F " << F ;
  std::cout << " rad " << pt.fit_rad;
  frags_at_popul = 1;
  print_timestamp();
  std::cout << std::endl;
}


MoverDE::MoverDE(boost::property_tree::ptree options, FitFunctionPtr scfxn_in, std::vector<Individual> initial_population) {
  srand(time(NULL));
  app_options = options;
  
  use_print_class = false;
  scfxn = scfxn_in;
  //  init_popul_ = init_popul_in;
#if(MPI_ENABLED)
  mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterRosettaCalculator(scfxn_in));

  std::string score_name = scfxn->name();
  std::string score_func = std::string(score_name.substr(score_name.rfind(" ") + 1));
  mpi_calculator->stage_name = score_func;
#endif


  mutation_strategy = std::string("default");
  boost::property_tree::ptree::const_assoc_iterator it = app_options.find("Protocol.mutation_strategy");
  if(app_options.not_found() == it) {
    mutation_strategy = app_options.get<std::string>("Protocol.mutation_strategy");
  }

  boost::optional< boost::property_tree::ptree& > child = app_options.get_child_optional( "Protocol.clean_gen_limit" );
  if( !child ) {
    clean_inds_gen_limit = 10000;
  } else {
    clean_inds_gen_limit =  app_options.get<int>("Protocol.clean_gen_limit");
  }
  child = app_options.get_child_optional( "Protocol.fragments_popul_module" );
  if( !child ) {
    local_search_popul_module = 1;
  } else {
    local_search_popul_module =  app_options.get<int>("Protocol.fragments_popul_module");
  }

  NP = app_options.get<int>("DE.NP");
  Gmax = 10000;
  D = scfxn->D();

  CR = app_options.get<double>("DE.CR");
  F = app_options.get<double>("DE.F");
  prot = app_options.get<std::string>("Protocol.prot");

  it = app_options.find("PrintPopul.gens");
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
  frags_at_popul = (app_options.get<bool>("Fragments.frags_at_popul")) ? 2 : 1 ;
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
  id = std::to_string( static_cast<int>(std::abs(rand())) );
  Parents parents;
  new_best_found = false;
  frags_at_popul = 2;
  std::cout << "init DE for function " << scfxn->name() << std::endl;
  popul = mpi_calculator->run(popul, 1);
  std::cout << "finish init popul " << std::endl;
  std::vector<Individual> trial_popul;
  popul_pdb.resize(NP);
  avg_acc = 0;
  for (int i = 0; i < popul.size(); i++) {
    avg_acc += popul[i].score;
    if (popul[i].score < popul[best_idx].score) {
      best_idx = i;
      best = popul[best_idx].score;
      if (best < previous_best) {
	previous_best = best;
      }
    }
  }

  std::cout << "at init " << best << " at " << best_idx << std::endl;

  
  while (gen_count < Gmax) {
    reset_stat();
    trial_popul.resize(0);
// #if(MPI_ENABLED)

//     popul = mpi_calculator->run(popul, frags_at_popul );
// #endif
    //calculate_distances_popul->build_pdb_population(popul,popul_pdb);
    for (int i = 0; i < NP; ++i) {
      select_parents(i, parents);
      Individual ind = sample_new_individual(i, parents, popul);
      trial_popul.push_back(ind);
      //trial_popul.push_back(popul[i]);
    }

#if(MPI_ENABLED)

    trial_popul = mpi_calculator->run(trial_popul, frags_at_popul );
#endif
    new_best_found = select_population(trial_popul);
    double c_best = (-1 * SCORE_ERROR_FIXED) + best;
    double c_avg = (-1 * SCORE_ERROR_FIXED) + (avg_acc / popul.size());
    print_generation_information(gen_count, new_best_found);
    ++gen_count;
    //    if ((std::abs(c_best) - std::abs(c_avg)) < 0.01)
    //  break;
  }


#if(USE_CRYO_EM)
  std::string prot_name = app_options.get<std::string>("Protocol.prot");
  std::string tag_dist_degrees = "deg" + app_options.get<std::string>("Protocol.disturb_degrees");
  std::string tag_reso = "res" + app_options.get<std::string>("Protocol.map_res");
  calculate_distances_popul->dump_popul(prot_name + "_" + tag_dist_degrees + "_" + tag_reso, popul);

#endif

}

void MoverDE::apply(std::vector<Individual>& popul_in, int NP_in, int Gmax_in) {
  gen_count = 0;
  int my_gen_count = 0;
  id = std::to_string( static_cast<int>(std::abs(rand())) );
  Parents parents;
  new_best_found = false;
  popul = popul_in;
  NP = NP_in;
  Gmax = Gmax_in;
  //std::cout << "INIT DE" << std::endl;
  double best_energy = popul_in[0].score;
  int best_idx_aux = 0;
  for (int i = 1; i < popul_in.size(); i++) {
    if (popul_in[i].score < best_energy) {
      best_idx_aux = i;
      best_energy = popul_in[i].score;
    }
  }
  best_idx = best_idx_aux;
  best = best_energy;

  popul_pdb.resize(NP);
  std::vector<Individual> trial_popul;
  while (my_gen_count < Gmax) {
    reset_stat();
    trial_popul.resize(0);
    calculate_distances_popul->build_pdb_population(popul_in,popul_pdb);
    for (int i = 0; i < NP; ++i) {
      select_parents(i, parents);
      Individual ind = sample_new_individual(i, parents, popul_in);
      trial_popul.push_back(ind);
    }

#if(MPI_ENABLED)
    trial_popul = mpi_calculator->run(trial_popul, frags_at_popul );
#endif
    new_best_found = select_population(trial_popul);
    popul_in = popul;
    // print_generation_information(gen_count, new_best_found);
    ++my_gen_count;
  }


#if(USE_CRYO_EM)
  std::string prot_name = app_options.get<std::string>("Protocol.prot");
  std::string tag_dist_degrees = "deg" + app_options.get<std::string>("Protocol.disturb_degrees");
  calculate_distances_popul->dump_popul(prot_name + "_" + tag_dist_degrees, popul);
#endif
}

void HybridMoverDE::apply_local_search_at_population() {
  if (gen_count % local_search_popul_module == 0) {
  local_search->reset_stats();
    double avg_before = 0, avg_after = 0;
    for (int i = 0; i < NP; ++i) {
      avg_before += (-1 * SCORE_ERROR_FIXED) + popul[i].score;
      local_search->apply(popul[i]);
      avg_after += (-1 * SCORE_ERROR_FIXED) + popul[i].score;
    }
    std::cout << local_search->print_stats() << " avg_before " << avg_before / NP << " avg_after " << avg_after / NP << std::endl;
  }
}


void SharedHybridMoverDE::apply_local_search_at_population() {
  if (gen_count % local_search_popul_module == 0) {
  local_search->reset_stats();
    double avg_before = 0, avg_after = 0;
    for (int i = 0; i < NP; ++i) {
      avg_before += (-1 * SCORE_ERROR_FIXED) + popul[i].score;
    }
    popul = mpi_calculator->run(popul);
    for (int i = 0; i < NP; ++i) {
      avg_after += (-1 * SCORE_ERROR_FIXED) + popul[i].score;
    }

    std::cout << local_search->print_stats() << " avg_before " << avg_before / NP << " avg_after " << avg_after / NP << std::endl;
  }
}

void HybridMoverDE::apply() {
  gen_count = 0;
  id = std::to_string( static_cast<int>(std::abs(rand())) );
  Parents parents;
  new_best_found = false;
  std::cout << "init Hybrid DE for function " << scfxn->name() << std::endl;
  std::vector<Individual> trial_popul;
  while (gen_count < Gmax) {
    reset_stat();
    scfxn->reset_stats();
    local_search->reset_stats();
    trial_popul.resize(0);
      //#pragma omp parallel for private(local_search)

    //apply_local_search_at_population();
    calculate_distances_popul->build_pdb_population(popul,popul_pdb);

    for (int i = 0; i < NP; ++i) {
      select_parents(i, parents);
      Individual ind = sample_new_individual(i, parents, popul);
      trial_popul.push_back(ind);
    }
    new_best_found = select_population(trial_popul);

    print_generation_information(gen_count, new_best_found);
    ++gen_count;
  }
}

bool MoverDE::select_population(const std::vector<Individual>& trial_popul) {
  // bool new_best_found = false;
  // std::vector<Individual> large_popul = trial_popul;
  // large_popul.insert(large_popul.end(), popul.begin(), popul.end());

  // std::sort(large_popul.begin(), large_popul.end(), popul_comparer);
  // large_popul.resize(popul.size());

  // if (large_popul[0].score < popul[0].score) {
  //   new_best_found = true;
  // }
  // popul = large_popul;
  // default version of selection
  // trial_sucess_n = 0;
  // avg_acc = 0;
  // for (int i = 0; i < NP; i++) {
  //   avg_acc += popul[i].score;
  //   if (trial_popul[i].score < popul[i].score) {
  //     popul[i] = trial_popul[i];
  //     trial_sucess_n++;
  //   }
  // }
  trial_sucess_n = 0;
  //  best_idx = 0;
  best = popul[best_idx].score;
  for (int i = 0; i < popul.size(); i++) {
    //std::cout << (-1 * SCORE_ERROR_FIXED) + popul[i].score << " ( "<< (-1 * SCORE_ERROR_FIXED) +  trial_popul[i].score << " ) ";
    if (trial_popul[i].score < popul[i].score) {

      popul[i] = trial_popul[i];
      if (popul[i].score < best) {
	best_idx = i;
	best = popul[i].score;
	trial_sucess_n++;
	if (best < previous_best) {
	  new_best_found = true;
	}

      }
    }
    avg_acc += popul[i].score;
  }
  //new_best_found = update_stat(0, trial_popul[0], popul) || new_best_found;
  return new_best_found;
}

bool CrowdingHybridMoverDE::select_population(const std::vector<Individual>& trial_popul) {

  return MoverDE::select_population(trial_popul);
  //  SharedHybridMoverDE::select_population(trial_popul);


  trial_sucess_n = 0;
  // Crowding DE: Rene Thomsen Algorithm
  // 1. Each trial can be replace only its nearest individual if the energy is better

  bool new_best_found = false;
  double previous_best = popul[best_idx].score;

  for (int i = 0; i < trial_popul.size(); i++) {
    int target_ind = 0;

    if (true) {
      target_ind = i;
    } else {
      target_ind = calculate_distances_popul->find_nearest(trial_popul[i], popul);
    }


    //std::cout << "select " << trial_popul[i].score << " < ? " << popul[target_ind].score << std::endl;
    if ( trial_popul[i].score < popul[target_ind].score ) {
      popul[target_ind] = trial_popul[i];
      new_best_found = true;
      trial_sucess_n++;
    }
  }

  std::cout << std::endl;
  avg_acc = 0;
  best_idx = 0;
  best = popul[0].score;
  //  std::cout << "check_popul ";
  for (int i = 0; i < popul.size(); i++) {
    //std::cout << (-1 * SCORE_ERROR_FIXED) + popul[i].score << " ( "<< popul[i].score << " ) ";
    avg_acc += popul[i].score;
    if (popul[i].score < popul[best_idx].score) {
      best_idx = i;
      best = popul[best_idx].score;
      if (best < previous_best) {
  	new_best_found = true;
      }
    }
  }
  //std::cout << std::endl;
  //std::cout << "best " << (-1 * SCORE_ERROR_FIXED) + best << " at " << best_idx << " " << popul[best_idx].score << std::endl;

  copy_popul_fitness_for_print = popul;

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
  return MoverDE::select_population(trial_popul);

  trial_sucess_n = 0;
  // Crowding DE: Rene Thomsen Algorithm

  bool new_best_found = false;
  double previous_best = popul[best_idx].score;
  double temp;
  std::vector<NeighStruct> aux(100);
  std::vector<double> aux_2(100);
  double value_dist_of_ind =  boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->distance_of_individual_partial_mario_euclidean(popul_pdb[best_idx],popul, temp,aux, aux_2);

  for (int i = 0; i < trial_popul.size(); i++) {
    int nearest_ind = calculate_distances_popul->find_nearest(trial_popul[i], popul);
    double previous_score = popul[nearest_ind].score;
    //std::cout << "nearest ind " << nearest_ind << " ind " << i << "(" << calculate_distances_popul->distance_between_inds(trial_popul[i], popul[nearest_ind])   << ") : " << trial_popul[i].score << " " << previous_score << std::endl;
    if ( trial_popul[i].score < popul[nearest_ind].score ) {
      popul[nearest_ind] = trial_popul[i];
      // int p;
      // std::cin >> p;
      new_best_found = true;
      trial_sucess_n++;
    }
  }

  avg_acc = 0;
  best_idx = 0;
  //std::cout << "check_popul " ;
  for (int i = 0; i < popul.size(); i++) {
    //std::cout <<  (-1 * SCORE_ERROR_FIXED) + popul[i].score << " ";
    avg_acc += popul[i].score;
    if (popul[i].score < popul[best_idx].score ) {
      best_idx = i;
      best = popul[best_idx].score;
      if (best < previous_best) {
	previous_best = best;
	new_best_found = true;
      }
    }
  }
  //std::cout << std::endl;
  //std::cout << "best " << best << " at " << best_idx << " " << popul[best_idx].score << std::endl;
  // std::string str_best_found = (new_best_found) ? "SI" : "NO";
  // std::cout << "best " << (-1 * SCORE_ERROR_FIXED) + best  << " idx " << best_idx << " previous " << previous_best << " " << str_best_found << std::endl;

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

  std::cout << "[GEN] " << gen_count << " " << (-1 * SCORE_ERROR_FIXED) + best << " " << (-1 * SCORE_ERROR_FIXED) + (avg_acc / popul.size()) << " at " << best_idx << " ";
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
  //  std::cout << local_search->print_stats() << std::endl;
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
#if(MPI_ENABLED)
  return 1000;
#else
  return scfxn->score(ind);
#endif
}

int
MoverDE::get_single_point_crossover(std::string ss) {
  int n_res = (rand() % static_cast<int>(ss.size()) )  + 1;
  while (ss[n_res] != 'L' || (n_res > (ss.size()) - 2) || (n_res < 2) ) {
    n_res = (rand() % static_cast<int>(ss.size())
 ) + 1;
  }
  int single_point_crossover = (n_res * 2) % D;
  return single_point_crossover;
}

Individual MoverDE::sample_new_individual(int i, const Parents& p, std::vector<Individual> popul) {
  Individual trial(size());
  int first_point = 0;
  int second_point = 0;
  int n_res = 0;
  int single_point_crossover = get_single_point_crossover(popul[i].ss);
  std::string ss = popul[0].ss;
  if (crossover_strategy == "loops") {
    std::vector<std::pair<int, int> > selected_regions = calculator_partial_rmsd->make_selected_regions();
    std::vector<std::pair<int, int> > loop_regions;
    for (int idx = 0; idx < selected_regions.size(); idx++) {
      if (ss[selected_regions[idx].first] == 'L') {
	loop_regions.push_back(selected_regions[idx]);
      }
    }
    std::pair<int, int> loop_cr = loop_regions[rand() % loop_regions.size() - 1];
    first_point = loop_cr.first * 2 - 1;
    second_point = loop_cr.second * 2 - 1;
  }

  for (int k = 0; k < D; ++k) {
    if (k > 1 && k % 2 != 0 ) n_res++;
    double jrand = URAND();
    bool mutation = true;
    if (crossover_strategy == "CR") {
      mutation = jrand < CR;
    } else {
      if (crossover_strategy == "SPC") {
	mutation = k < single_point_crossover;
      }
    }
    bool mutation_only_loops = false;
    if (mutation_strategy == "only_loops") {
      mutation_only_loops = true;
    }

    if (crossover_strategy == "loops") {
      if ( (k >= first_point) && (k <= second_point)) {
	trial.vars[k] = popul[i].vars[k];
	mutation = false;
      }
    }

    if ( mutation ) {
      double diff =  popul[p.x2].vars[k] - popul[p.x3].vars[k];
      // diff_avg_bef += std::abs(diff);
      // while (diff > 1) { diff -= 1; }
      // while (diff < -1) { diff += 1; }
      // diff_avg_aft += std::abs(diff);
      // if (mutation_only_loops) {
      // 	if  (popul[p.x1].ss[n_res] == 'L')   {
      // 	  trial.vars[k] = popul[p.x1].vars[k] + F * (diff);
      // 	} else {
      // 	  trial.vars[k] = popul[p.x1].vars[k] + 0 * (diff);
      // 	}
      // } else {
      trial.vars[k] = popul[p.x1].vars[k] + F * (diff);

      //	trial.vars[k] = popul[i].vars[k];
	//}

      // //           std::cout << "n_res " << n_res << " " << popul[i].ss[n_res] <<  std::endl;
      // if (trial.vars[k] > DE_LIMIT) {
      // 	trial.vars[k] = DE_LIMIT - std::abs(trial.vars[k] - DE_LIMIT);
      // }
      // if (trial.vars[k] < (-1*DE_LIMIT)) {
      // 	trial.vars[k] = (-1) + std::abs(trial.vars[k] - (-1*DE_LIMIT));
      // }

    } else {
      if (crossover_strategy != "loops") {
	trial.vars[k] = popul[i].vars[k];
      }
    }
  }


  //trial.vars = popul[p.x1].vars;
  if (CR > 0.5) {
    trial.omega = popul[p.x1].omega;
    trial.ss = popul[p.x1].ss;
  } else {
    trial.omega = popul[i].omega;
    trial.ss = popul[i].ss;
  }
  double result = score(trial);
  trial.score = popul[i].score;
  return trial;
}

int MoverDE::size() { return D; }

int MoverDE::get_best_rand_ind(int target, std::string select_parents_strategy) {
   int module_val = 5;
    if (select_parents_strategy == "best_ten") {
      module_val = 10;
    }
    std::vector<std::pair<double, int> > a;

    for (int i = 0 ; i < NP ; i++) {
      a.push_back (make_pair (popul[i].score,i)); // k = value, i = original index
    }
    std::sort(a.begin(),a.end());
    int rand_five = rand() % module_val;
    do {
      rand_five = rand() % module_val;
    } while ( a[rand_five].second == target );

    return rand_five;
}

void MoverDE::select_parents(int i, Parents& parent) {
  if (select_parents_strategy == "best") {
    parent.x1 = best_idx;
  }
  if (select_parents_strategy == "rand") {
    if (URAND() < 0.5) {
      parent.x1 = best_idx;
    } else {
      do {
	parent.x1 = rand() % NP;
      } while (parent.x1 == i);
    }
  }
  if ( (select_parents_strategy == "best_five") || ( select_parents_strategy == "best_ten") )  {
    int best_rand_ind = get_best_rand_ind(i, select_parents_strategy );
    parent.x1 = best_rand_ind;
  }
  if (select_parents_strategy == "nearest")  {
    int nearest_ind = calculate_distances_popul->find_nearest(popul[i], popul, i, popul_pdb);
    parent.x1 = nearest_ind;
    if (parent.x1 == i) {
      std::cout << "ERROR EN NEAREST" << std::endl;
      int p;
      std::cin >> p;
    }
  }

  if (select_parents_strategy == "top_nearest") {
   double temp;
   std::vector<NeighStruct> aux(100);
   std::vector<double> aux_2(100);
   double value_dist_of_ind =  boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->distance_of_individual_partial_mario_euclidean(popul_pdb[0],popul, temp,aux, aux_2);
    std::vector<int> nearest_top_ind = calculate_distances_popul->find_top_nearest(popul[i], popul, i, popul_pdb);

    int tam = nearest_top_ind.size();
    do {
      parent.x1 = nearest_top_ind[ rand() % tam ];
    } while (parent.x1 == i);
    do {
      parent.x2 = nearest_top_ind[ rand() % tam ];
    } while ((parent.x2 ==i) || (parent.x2 == parent.x1));

    do {
      parent.x3 = nearest_top_ind[ rand() % tam ];
    } while ((parent.x3 == i) || (parent.x3 == parent.x2) || (parent.x3 == parent.x1));


   //    std::cout << "x1 " << parent.x1 << " x2 " << parent.x2 << " x3 " << parent.x3 << std::endl;

  } else {
    do {
      parent.x2 = rand() % NP;
    } while ((parent.x2 ==i) || (parent.x2 == parent.x1));

    do {
      parent.x3 = rand() % NP;
    } while ((parent.x3 == i) || (parent.x3 == parent.x2) || (parent.x3 == parent.x1));

  }

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
    // int p;
    // std::cin >> p;
  }
}

void CrowdingHybridMoverDE::print_fitness_population(int gen_count) {
  //if ((gen_count < 25) || (gen_count % 50 == 0)) {
  if (true) {
    std::cout << "[POP] ";
    for (int i = 0; i < copy_popul_fitness_for_print.size(); i++) {
      std::cout << (-1 * SCORE_ERROR_FIXED) + copy_popul_fitness_for_print[i].score << " ";
    }
    std::cout << std::endl;

    // for (int i = 0; i < copy_popul_fitness_for_print.size(); i++) {
    //   std::cout << "ind " << i << " : " << copy_popul_fitness_for_print[i].gen_acc << std::endl;
    // }

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

  // std::cout << "[DIST_POP] ";
  // std::cout << ind_distance_avg[best_idx]  + (-1 * SCORE_ERROR_FIXED) << " ";
  // for (int i = 0; i < large_NP; i++) {
  //   std::cout << ind_distance_avg[i] + (-1 * SCORE_ERROR_FIXED) << " ";
  // }
  // std::cout << std::endl;

  // InterDistancesGraph inter_distances_calc;
  // /* inter distances graphic */
  // double min_inter_distance = 1000, max_inter_distance = -1000;
  // for (int i = 0; i < ind_inter_distance.size(); i++) {
  //   // if (min_inter_distance > ind_inter_distance[i]) {
  //   //   min_inter_distance = ind_inter_distance[i];
  //   // }
  //   // if (max_inter_distance < ind_inter_distance[i]) {
  //   //   max_inter_distance =ind_inter_distance[i];
  //   // }
  //   //    std::cout << ind_inter_distance[i] << std::endl;
  //   inter_distances_calc.add_occurrence(ind_inter_distance[i]);
  // }
  // //  std::cout << "min_inter_distance " << min_inter_distance << " max_inter_distance " << max_inter_distance << std::endl;
  // //std::cout << "[INTER] ";
  // //std::cout << inter_distances_calc.print_graph() << std::endl;
  // /* END inter distances graphic */

  // if (neighs_per_ind.size() > 0) {
  //   std::cout << "[NEIGH] ";

  //   for (int i = 0; i < large_NP; i++) {
  //     std::cout << neighs_per_ind[i].size() << " ";
  //   }
  //   std::cout << std::endl;
  // } else {
  //   std::cout << "[NEIGH] ";
  //   for (int i = 0; i < NP; i++) {
  //     std::cout << result.neigh_per_ind[i].size() << " ";
  //   }
  //   std::cout << std::endl;
  // }

}


void MoverDE::print_distances_population() {
  // CalculateDistancePopulation::DistancesResult result =
  //                                       calculate_distances_popul->run(popul);

  // std::vector<double> ind_distance_avg = result.distances_of_population;
 // std::vector<double> rmsd_native = result.rmsd_to_native;
  std::vector<double> rmsd_native = calculate_distances_popul->run_rmsd_to_native(popul);
  // std::vector<double> ind_inter_distance = result.distances_of_population;

  std::cout << "[RMSD_NATIVE] ";
  for (int i = 0; i < popul.size(); i++) {
    std::cout << rmsd_native[i] << " ";
  }
  std::cout << std::endl;

  // std::cout << "[DIST_POP] ";
  // std::cout << " ";
  // //for (int i = 0; i < ind_distance_avg.size(); i++) {
  // for (int i = 0; i < 50; i++) {
  //   std::cout << "1 ";
  // }
  // std::cout << std::endl;

  // InterDistancesGraph inter_distances_calc;

  // //for (int i = 0; i < ind_inter_distance.size(); i++) {
  //   for (int i = 0; i < 10; i++) {
  //           inter_distances_calc.add_occurrence(0);
  //   }

  //   std::cout << "[INTER] ";
  //   std::cout << inter_distances_calc.print_graph() << std::endl;

  std::vector<std::vector<double> > rmsd_native_all = calculator_partial_rmsd->partial_rmsd_of_population(popul);
  std::string head = calculator_partial_rmsd->make_head(popul);
  // [ {zona1}, {zona2}, {zona3}, {zona4} ],
  // [ rmsd_native_all[0][zona_i].... ],
  // [ rmsd_native_all[NP][zona_i].... ],
  /*  std::cout << "[PARTIAL_RMSD] [ ";
  std::cout << head << " , ";
  int n_pieces = rmsd_native_all[0].size();
  for (int i = 0; i < rmsd_native_all.size(); i++) {
    std::cout << " [ ";
    for (int j = 0; j < n_pieces; j++) {
      if (j != n_pieces - 1) {
	std::cout << rmsd_native_all[i][j] << " , ";
      } else {
	std::cout << rmsd_native_all[i][j] << " ] ";
      }
    }
    if (i != rmsd_native_all.size() - 1) {
      std::cout << " , ";
    }
  }
  std::cout << " ] " << std::endl;
  */
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
    //if ((new_best_found) || (gen_count == 0)) {
    if ( (gen_count % 10) == 0  ) {
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


      std::string url_best = print_best_ind->print(popul[ind_index], id);
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
  //  avg_acc += popul[i].score;
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
  id = std::to_string( static_cast<int>(std::abs(rand())) );
  Parents parents;
  new_best_found = false;
  std::cout << "init Shared Hybrid DE for function " << scfxn->name() << std::endl;
  std::vector<Individual> trial_popul;
  double global_min_rmsd, global_max_rmsd, global_avg_rmsd;
  while (gen_count < Gmax) {
    timestamp_t t_ini = get_timestamp();
    reset_stat();
    local_search->reset_stats();
    trial_popul.resize(0);
    global_min_rmsd = 100000;
    global_max_rmsd = -100000;
    global_avg_rmsd = 0;

    double total_diff_avg_bef=0, total_diff_avg_aft = 0;
    diff_avg_bef=0; diff_avg_aft = 0;

    int cnt = 0;
    double duration_acc = 0;

    apply_local_search_at_population();
    calculate_distances_popul->build_pdb_population(popul,popul_pdb);
    for (int i = 0; i < NP; ++i) {
      select_parents(i, parents);
      Individual ind = sample_new_individual(i, parents, popul);
      trial_popul.push_back(ind);
    }

    timestamp_t t0 = get_timestamp();


    trial_popul = mpi_calculator->run(trial_popul, frags_at_popul );
    timestamp_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;
    std::cout << "evaluation stage time" << secs << std::endl;

    new_best_found = select_population(trial_popul);

    print_generation_information(gen_count, new_best_found);
    ++gen_count;

    timestamp_t t_final = get_timestamp();

   double secs_t = (t_ini - t_final) / 1000000.0L;
    std::cout << "gen stage time" << secs_t << std::endl;
  }
}


void ResetOldIndsHybridDE::apply() {
  gen_count = 0;
  id = std::to_string( static_cast<int>(std::abs(rand())) );
  Parents parents;
  new_best_found = false;
  std::cout << "init Reset Old DE for function " << scfxn->name() << std::endl;
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

    int cnt = 0;
    double duration_acc = 0;
    old_population = popul;

    apply_local_search_at_population();
    calculate_distances_popul->build_pdb_population(popul,popul_pdb);
    for (int i = 0; i < NP; ++i) {
      select_parents(i, parents);
      Individual ind = sample_new_individual(i, parents, popul);
      trial_popul.push_back(ind);
    }

    new_best_found = select_population(trial_popul);

    clean_old_inds();

    print_generation_information(gen_count, new_best_found);

    ++gen_count;
  }
}

void ResetOldIndsHybridDE::clean_old_inds() {
    // clean old individuals
  if (clean_inds_gen_limit < 100) {
    for (int i = 0; i < NP; i++) {
      if ( (i != best_idx) && ( popul[i].gen_acc >= clean_inds_gen_limit ) ) {
	//if ( (i != best_idx) && ( popul[i].gen_acc >= clean_inds_gen_limit ) ) {
	double prev_score = popul[i].score;
	greedy_mover->apply(popul[i]);
	popul[i].gen_acc = 0;
	//	std::cout << "!!! RESET IND " << i  << " : " << (-1*SCORE_ERROR_FIXED) + prev_score << " to " <<  (-1*SCORE_ERROR_FIXED) +  popul[i].score << std::endl;
      }
    }
  }
}

bool
ResetOldIndsHybridDE::select_population(const std::vector<Individual>& trial_popul) {
  trial_sucess_n = 0;
 for (int i = 0; i < NP; i++) {
    if (trial_popul[i].score < popul[i].score) {
      popul[i] = trial_popul[i];
      popul[i].gen_acc = 0;
      trial_sucess_n++;
      new_best_found = true;
    }
  }

  best_idx = 0;
  best = popul[best_idx].score;
  for (int i = 0; i < popul.size(); i++) {
    //std::cout << (-1 * SCORE_ERROR_FIXED) + popul[i].score << " ( "<< popul[i].score << " ) ";
    avg_acc += popul[i].score;
    if (popul[i].score < popul[best_idx].score) {
      best_idx = i;
      best = popul[best_idx].score;
      if (best < previous_best) {
	previous_best = best;
	new_best_found = true;
      }
    }
  }

  return new_best_found;
}

bool
ResetOldIndsCrowdingHybridDE::select_population(const std::vector<Individual>& trial_popul) {
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
      popul[nearest_ind].gen_acc = 0;
    }
  }
  for (int i = 0; i < NP; i++) {
    if (popul[i].score == old_population[i].score) {
      popul[i].gen_acc = popul[i].gen_acc + 1;
    }
  }
  avg_acc = 0;
  //  best_idx = 0;
  best = popul[best_idx].score;
  //  std::cout << "check_popul ";
  for (int i = 0; i < popul.size(); i++) {
    //std::cout << (-1 * SCORE_ERROR_FIXED) + popul[i].score << " ( "<< popul[i].score << " ) ";
    avg_acc += popul[i].score;
    if (popul[i].score < popul[best_idx].score) {
      best_idx = i;
      best = popul[best_idx].score;
      if (best < previous_best) {
	previous_best = best;
	new_best_found = true;
      }
    }
  }
  //std::cout << std::endl;
  //std::cout << "best " << (-1 * SCORE_ERROR_FIXED) + best << " at " << best_idx << " " << popul[best_idx].score << std::endl;
  return new_best_found;
}

void CrowdingHybridMoverDE::print_distances_population() {
  MoverDE::print_distances_population();
}

// otra forma para sample individual
/*      if (trial.vars[k] > DE_LIMIT) {
	trial.vars[k] = 1.0 - (URAND() * ( 1.0 * 0.05));
      }
      if (trial.vars[k] < (-1*DE_LIMIT)) {
	trial.vars[k] = -1.0 + (URAND() * ( -1.0 * 0.05));
      }


 */

#if(MPI_ENABLED)
void MPICrowdingMoverDE::apply_local_search_at_population() {
  if ((local_search_popul_module != 0) && ( gen_count % local_search_popul_module == 0) ) {
    local_search->reset_stats();
    double avg_before = 0;
    for (int i = 0; i < NP; ++i) {
      avg_before += (-1 * SCORE_ERROR_FIXED) + popul[i].score;
    }
    popul = mpi_calculator->run(popul);
    double avg_after= 0;
    for (int i = 0; i < NP; ++i) {
      avg_after += (-1 * SCORE_ERROR_FIXED) + popul[i].score;
    }
    std::cout << local_search->print_stats() << " avg_before " << avg_before / NP << " avg_after " << avg_after / NP << std::endl;
  }
}

void MPICrowdingMoverDE::analysis_of_population(std::vector<Individual> popul, std::vector<core::pose::PoseOP> popul_pdb){

  int best_ind = 0;
  int worst_ind = 0;
  double best_value =popul[best_ind].score;
  double worst_value = popul[worst_ind].score;
  for (int i = 0; i < popul.size(); i++) {
    if (popul[i].score < best_value) {
      best_value = popul[i].score;
      best_ind = i;
    }
     if (popul[i].score > worst_value) {
      worst_value = popul[i].score;
      worst_ind = i;
    }
  }

   boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->build_pdb_population(popul, popul_pdb);
   boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->build_inter_distances_of_population(popul);
   double temp;
   std::vector<NeighStruct> aux(100);
   std::vector<double> aux_2(100);
   double value_dist_of_ind =  boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->distance_of_individual_partial_mario_euclidean(popul_pdb[best_ind],popul, temp,aux, aux_2);
   //std::cout << "value dist " << value_dist_of_ind << std::endl;
   //pose_distance_calculation(core::pose::Pose& pose1, core::pose::Pose& pose2)

   double best_vs_worst_distance = boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->pose_distance_calculation(popul_pdb[best_ind],popul_pdb[worst_ind]);

   double best_native_distance = boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->pose_distance_calculation(popul_pdb[best_ind],calculate_distances_popul->native_);

   double worst_native_distance = boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->pose_distance_calculation(calculate_distances_popul->native_,popul_pdb[worst_ind]);

   std::cout << "[DIST_POPUL] ";
   std::cout << "best ind " << best_ind << " " <<  (-1 * SCORE_ERROR_FIXED) + best_value << " ";
   std::cout << "worst ind " << worst_ind << " " <<  (-1 * SCORE_ERROR_FIXED) + worst_value << " ";
   std::cout << "best to native " << best_native_distance << " ";
   std::cout << "worst to native " << worst_native_distance << " ";
   std::cout << "best to worst " << best_vs_worst_distance << " ";
   /*double diff =  std::abs(best_native_distance - worst_native_distance);

   if (diff > 0) {
       current_cde_radius = diff;
   } else {
	current_cde_radius = 1;
   }*/

   //current_cde_radius = 0.00303504;
   //   current_cde_radius = 0.03;
   std::cout << std::endl;
   //   exit(1);
}


std::vector<std::vector<Individual> >
MPISeedsMoverDE::create_seeds(std::vector<Individual> popul, std::vector<core::pose::PoseOP> popul_pdb) {
  std::vector<IndIndex> sort_popul_index;
  for (int i = 0; i < popul.size(); i++) {
    sort_popul_index.push_back(IndIndex(i, popul[i]));
  }

  std::vector<Individual> sort_popul = popul;
  std::sort(sort_popul.begin(), sort_popul.end(), popul_comparer );
  PopulSeedsComparer popul_seeds_comparer;
  std::sort(sort_popul_index.begin(), sort_popul_index.end(), popul_seeds_comparer );

  calculate_distances_popul->build_pdb_population(sort_popul, popul_pdb);
  if (app_options.get<std::string>("Protocol.distance_strategy") == "euclidean_partial_mario") {
  boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->build_inter_distances_of_population(sort_popul);
  }

  std::vector<int> current_seed;
  current_seed.push_back(sort_popul_index[0].index_);

  //current_cde_radius = 0.03;
  current_cde_radius = 7.0;
  //  current_cde_radius = 0.01;
  std::map<int, std::vector<int> > population_per_seed;
  std::vector<int> processed_inds(sort_popul.size(), 0);
  population_per_seed[current_seed[0]] = std::vector<int>();
  population_per_seed[current_seed[0]].push_back( current_seed[0] );
  for (int i = 1; i < sort_popul.size(); i++) {
    if (processed_inds[sort_popul_index[i].index_] == 0) {
      bool found = false;
      for (int j = 0; j < current_seed.size() && !found; j++) {
	double dist = 0;
	if (app_options.get<std::string>("Protocol.distance_strategy") == "euclidean_partial_mario") {
	  dist = boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->pose_distance_calculation(popul_pdb[current_seed[j]],popul_pdb[i]);
	} else {
	  dist = calculate_distances_popul->current_distance_calculation(*popul_pdb[current_seed[j]],*popul_pdb[i]);
	}
	if ( dist <= current_cde_radius) {
	  found = true;
	  population_per_seed[current_seed[j]].push_back(sort_popul_index[i].index_);
	  processed_inds[sort_popul_index[i].index_]  = 1;
	  break;
	}
      }
      if (!found) {
	current_seed.push_back(sort_popul_index[i].index_);
	processed_inds[sort_popul_index[i].index_]  = 1;
	population_per_seed[sort_popul_index[i].index_] = std::vector<int>();
	population_per_seed[sort_popul_index[i].index_].push_back(sort_popul_index[i].index_);
      }
    }
  }

  int min = 1000;
  int min_index = 0;
  // Limit current seeds to 3
  int max_subpopuls = 3;
  if (current_seed.size() > max_subpopuls) {
    for (int i = max_subpopuls; i < current_seed.size(); i++) {
      for (int j = 0; j < population_per_seed[current_seed[i]].size(); j++) {
	for (int k = 0; k < max_subpopuls; k++) {
	  double dist = 0;
	  if (app_options.get<std::string>("Protocol.distance_strategy") == "euclidean_partial_mario") {
	    dist = boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->pose_distance_calculation(popul_pdb[current_seed[k]], popul_pdb[population_per_seed[current_seed[i]][j]] );
	  } else {
	    dist = calculate_distances_popul->current_distance_calculation(*popul_pdb[current_seed[k]], *popul_pdb[population_per_seed[current_seed[i]][j]] );
	  }

	  if (dist < min) {
	    min = dist;
	    min_index = k;
	  }
	}
	population_per_seed[current_seed[min_index]].push_back(population_per_seed[current_seed[i]][j]);
      }
    }
    current_seed.resize(max_subpopuls);
  }


  //std::cout << std::endl;
  std::cout << "[CURRENT_SEEDS] ";
  for (int i = 0; i < current_seed.size(); i++) {
    std::cout << i << " : " << popul[current_seed[i]].score << " , ";
  }
  std::cout << std::endl;

  std::vector<std::vector<Individual> > subpopuls;
  int count = 0;
  for (int i = 0; i < current_seed.size(); i++) {
    //     std::cout << "SEED " << current_seed[i] << " : " <<  popul[current_seed[i]].score << std::endl;

     subpopuls.push_back(std::vector<Individual>());
    for (int j = 0; j < population_per_seed[current_seed[i]].size(); j++) {
      //std::cout << "\t" << population_per_seed[current_seed[i]][j] << " : " <<  popul[population_per_seed[current_seed[i]][j]].score << std::endl ;
      subpopuls[i].push_back(popul[population_per_seed[current_seed[i]][j]]);
      count++;
    }
  }
  //std::cout << std::endl;
  //  std::cout << "total processed elements " << count << std::endl;
  for (int i = 0; i < subpopuls.size(); i++) {

    //subpopuls[i] = mpi_calculator->run( subpopuls[i], 1 );
    std::vector<double> rmsd_subpopul = calculate_distances_popul->run_rmsd_to_native(subpopuls[i]);
    std::cout << "subpopul " << i << " energy [";
    for (int j = 0; j < subpopuls[i].size(); j++) {
      if (j == subpopuls[i].size() - 1) {
	std::cout << " " << subpopuls[i][j].score << " ";
      }else {
	std::cout << " " << subpopuls[i][j].score << " ,";
      }
    }

    std::cout << "]" << std::endl;

    std::cout << "subpopul " << i << " rmsd [";
    for (int j = 0; j < subpopuls[i].size(); j++) {
      if (j == subpopuls[i].size() - 1) {
	std::cout << " " << rmsd_subpopul[j] << " ";
      }else {
	std::cout << " " << rmsd_subpopul[j] << " ,";
      }
    }
    std::cout << "]" << std::endl;
  }

  //create subpopuls
  int m_size = popul.size() * 2;
  int ind_size = popul[0].vars.size();
  for (int i = 0; i < subpopuls.size(); i++) {
    while (subpopuls[i].size() < m_size) {
      Individual nuevo(ind_size);
      init_popul_->disturb_individual(nuevo, ind_size);
      subpopuls[i].push_back(nuevo);
    }

    for (int j = 1; j < subpopuls[i].size(); j++) {
      if (std::abs(subpopuls[i][0].score - subpopuls[i][j].score) > 0.01) {
	 init_popul_->disturb_individual(subpopuls[i][j], ind_size);
      }
    }
  }

  // score & print subpopuls
  for (int i = 0; i < subpopuls.size(); i++) {
    subpopuls[i] = mpi_calculator->run( subpopuls[i], 1 );
  }

  return subpopuls;
}


void MPISeedsMoverDE::apply() {
  gen_count = 0;
  id = std::to_string( static_cast<int>(std::abs(rand())) );
  Parents parents;
  new_best_found = false;
  std::cout << "init MPI Seeds Mover " << scfxn->name() << std::endl;
  std::vector<Individual> trial_popul;
  std::vector<Individual> main_popul = popul;
  int Gmax_in = 20;
  int my_gen_count = 0;
  int main_NP = NP;
  std::vector<Individual> result_popul;
  while (my_gen_count < Gmax_in) {
    result_popul.resize(0);
    reset_stat();
    trial_popul.resize(0);
    apply_local_search_at_population();
    calculate_distances_popul->build_pdb_population(popul,popul_pdb);
    //analysis_of_population(popul, popul_pdb);
    std::vector<std::vector<Individual> > subpopuls = create_seeds(popul, popul_pdb);

    std::vector<Individual> super_popul;
    result_popul.resize(0);
    for (int i = 0; i < subpopuls.size(); i++) {
      NP = main_NP;
      //std::cout << "start for subpopul " << i << " size " << subpopuls[i].size() << " best " << subpopuls[i][0].score << std::endl;
      MoverDE::apply(subpopuls[i], subpopuls[i].size(), 5 );
      std::sort(subpopuls[i].begin(), subpopuls[i].end(), popul_comparer);
      //std::cout << "subpopul " << i << " : " << subpopuls[i][0].score << " , " << subpopuls[i][subpopuls[i].size() - 1].score << std::endl;
      // result_popul.push_back(subpopuls[i][0]);
      // subpopuls[i].erase(subpopuls[i].begin());
      for (int j = 0; j < subpopuls[i].size(); j++) {
	super_popul.push_back(subpopuls[i][j]);
      }
    }

    std::sort(super_popul.begin(), super_popul.end(), popul_comparer);
    // // for (int i = result_popul.size(); i < main_NP; i++) {
    // //   result_popul.push_back(super_popul[i]);
    // // }
    result_popul = super_popul;
    // // takes one ind per subpopul until completes the main popul
    // while(result_popul.size() < main_NP) {
    //   for (int i = 0; i < subpopuls.size() && result_popul.size() < main_NP; i++) {
    // 	result_popul.push_back(subpopuls[i][0]);
    // 	subpopuls[i].erase(subpopuls[i].begin());
    //   }
    // }
    // std::sort(result_popul.begin(), result_popul.end(), popul_comparer);
    best_idx = 0;
    if (popul[best_idx].score > result_popul[0].score) {
      new_best_found = true;
    }
    result_popul.resize(main_NP);
    main_popul = result_popul;
    popul = result_popul;
    copy_popul_fitness_for_print = popul;
    best = popul[best_idx].score;
    avg_acc  = 0;
    for (int i = 0; i < popul.size(); i++) {
      avg_acc += popul[i].score;
    }

    print_generation_information(my_gen_count, new_best_found);
    ++my_gen_count;
    ++gen_count;
  }
}

void MPICrowdingMoverDE::apply() {
  gen_count = 0;
  id = std::to_string( static_cast<int>(std::abs(rand())) );
  Parents parents;
  new_best_found = false;
  std::cout << "init MPI Crowding Mover " << scfxn->name() << std::endl;
  std::vector<Individual> trial_popul;
  double global_min_rmsd, global_max_rmsd, global_avg_rmsd;
  bool equal_value = false;

  while (gen_count < Gmax && !equal_value) {
    timestamp_t t_ini = get_timestamp();
    reset_stat();
    trial_popul.resize(0);
    global_min_rmsd = 100000;
    global_max_rmsd = -100000;
    global_avg_rmsd = 0;
    //apply_local_search_at_population();

    calculate_distances_popul->build_pdb_population(popul,popul_pdb);
    //analysis_of_population(popul, popul_pdb);

    timestamp_t t_sample = get_timestamp();
    for (int i = 0; i < NP; ++i) {
      select_parents(i, parents);
      Individual ind = sample_new_individual(i, parents, popul);
      trial_popul.push_back(ind);
    }
    timestamp_t t_sample_end = get_timestamp();
    double secs_sample = (t_sample_end - t_sample) / 1000000.0L;
    std::cout << "sample stage time " << secs_sample << std::endl;

    timestamp_t t0 = get_timestamp();


    trial_popul = mpi_calculator->run(trial_popul, frags_at_popul );
    timestamp_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;
    std::cout << "evaluation stage time " << secs << std::endl;


    timestamp_t t_start_select = get_timestamp();
    new_best_found = select_population(trial_popul);
    timestamp_t t_end_select = get_timestamp();
    double secs_selection = (t_end_select - t_start_select) / 1000000.0L;
    std::cout << "selection stage time " << secs_selection << std::endl;

    double avg_acc = 0;
    for (int p = 0; p < popul.size(); p++) {
      avg_acc += popul[p].score;
    }
    avg_acc = avg_acc / popul.size();

    if (std::abs(avg_acc - popul[best_idx].score) < 0.1) {
      //equal_value = true;
    }

    if (false) {
      // quiero imprimir el score vs rmsd de cada Individuo que me importa para mostrarlos en una grafica.
      double best_rmsd = 0;
      int nearest = calculate_distances_popul->find_nearest(trial_popul[best_idx], popul);
      std::vector<double> rmsd_native_popul = calculate_distances_popul->run_rmsd_to_native(popul);
      std::vector<double> rmsd_native_trial = calculate_distances_popul->run_rmsd_to_native(trial_popul);
      std::cout << "[BEST_IND_ANALYSIS] { best : ( " << (-1 * SCORE_ERROR_FIXED) + popul[best_idx].score << ", " << rmsd_native_popul[best_idx] << " ) , " << " ";
      std::cout << "trial : (" << (-1 * SCORE_ERROR_FIXED) + trial_popul[best_idx].score << " , " << rmsd_native_trial[best_idx] << ") , " << " ";
      std::cout << "parents : [ (" << (-1 * SCORE_ERROR_FIXED) + popul[parents.x1].score <<"," <<rmsd_native_popul[parents.x1] << " ) , "  << "(" << (-1 * SCORE_ERROR_FIXED) + popul[parents.x2].score <<"," <<rmsd_native_popul[parents.x2] << " ), " << "(" << (-1 * SCORE_ERROR_FIXED) + popul[parents.x3].score <<"," <<rmsd_native_popul[parents.x3] << " ) ] , "  << " ";
      std::cout << "nearest : (" << (-1 * SCORE_ERROR_FIXED) + popul[nearest].score << "," << rmsd_native_popul[nearest]  << " ) } " << std::endl;

    }
    print_generation_information(gen_count, new_best_found);
    ++gen_count;

    timestamp_t t_final = get_timestamp();
    double secs_t = (t_final - t_ini) / 1000000.0L;
    std::cout << "gen stage time " << secs_t << std::endl;

  }
}

double SharedHybridMoverDE::score(Individual& ind) {
  return 1;
}



double MPICrowdingMoverDE::score(Individual& ind) {
  return 1;
}

double MPIResetOldCrowdingHybridDE::score(Individual& ind) {
  return 1;
}

void MPIResetOldCrowdingHybridDE::apply() {
  gen_count = 0;
  id = std::to_string( static_cast<int>(std::abs(rand())) );
  Parents parents;
  new_best_found = false;
  std::cout << "init MPI Reset Old Crowding Mover " << scfxn->name() << std::endl;
  std::vector<Individual> trial_popul;
  double global_min_rmsd, global_max_rmsd, global_avg_rmsd;
  while (gen_count < Gmax) {
    reset_stat();
    trial_popul.resize(0);
    global_min_rmsd = 100000;
    global_max_rmsd = -100000;
    global_avg_rmsd = 0;

    old_population = popul;
    apply_local_search_at_population();
    for (int i = 0; i < NP; ++i) {
      if (popul[i].score < old_population[i].score ) {
	popul[i].gen_acc = 0;
      }
    }
    old_population = popul;

    calculate_distances_popul->build_pdb_population(popul,popul_pdb);
    for (int i = 0; i < NP; ++i) {
      select_parents(i, parents);
      Individual ind = sample_new_individual(i, parents, popul);
      trial_popul.push_back(ind);
    }

    trial_popul = mpi_calculator->run(trial_popul);
    new_best_found = select_population(trial_popul);

    clean_old_inds();

    print_generation_information(gen_count, new_best_found);
    ++gen_count;
  }
}




void MPIfastCrowdingMoverDE::apply() {
  gen_count = 0;
  id = std::to_string( static_cast<int>(std::abs(rand())) );
  Parents parents;
  new_best_found = false;
  std::cout << "init MPI super fast Crowding Mover " << scfxn->name() << std::endl;
  std::vector<Individual> trial_popul;
  double global_min_rmsd, global_max_rmsd, global_avg_rmsd;
  bool equal_value = false;
  //mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new  MasterScatterGather(scfxn, calculate_distances_popul));
  mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterEvaluateAndNearestCalculator(scfxn, calculate_distances_popul));

  while (gen_count < Gmax && !equal_value) {
    timestamp_t t_ini = get_timestamp();
    reset_stat();
    trial_popul.resize(0);
    global_min_rmsd = 100000;
    global_max_rmsd = -100000;
    global_avg_rmsd = 0;
    calculate_distances_popul->build_pdb_population(popul,popul_pdb);
    timestamp_t t_sample = get_timestamp();
    for (int i = 0; i < NP; ++i) {
      select_parents(i, parents);
      Individual ind = sample_new_individual(i, parents, popul);
      trial_popul.push_back(ind);
    }
    timestamp_t t_sample_end = get_timestamp();
    double secs_sample = (t_sample_end - t_sample) / 1000000.0L;
    std::cout << "sample stage time " << secs_sample << std::endl;

    timestamp_t t0 = get_timestamp();

    bool test_gecco2015_only_score_individuals = false;
    std::vector<IndMPI> result_from_mpi;
    if (test_gecco2015_only_score_individuals) {
      result_from_mpi = mpi_calculator->run(popul, trial_popul, 1);
    } else {
      result_from_mpi = mpi_calculator->run(popul, trial_popul, 2);
    }
    nearest_inds.resize(NP);
    for (int k= 0; k< result_from_mpi.size(); k++) {
      trial_popul[k] = result_from_mpi[k].ind;
      nearest_inds[k] = result_from_mpi[k].nearest;
    }

    std::vector<int> check_nearest;
    for (int k = 0; k < trial_popul.size(); k++) {
	int target_ind = calculate_distances_popul->find_nearest(trial_popul[k], popul);
	check_nearest.push_back(target_ind);
    }

    //  for (int k = 0; k < trial_popul.size(); k++) {
    //    if (check_nearest[k] != nearest_inds[k]) {
    // 	 std::cout << "DIFFERENT NEAREST " << k << " : " << check_nearest[k] << " != " << nearest_inds[k] << std::endl;
    //    }
    // }

    timestamp_t t1 = get_timestamp();
    double secs = (t1 - t0) / 1000000.0L;
    std::cout << "evaluation stage time " << secs << std::endl;


    timestamp_t t_start_select = get_timestamp();
    new_best_found = select_population(trial_popul);
    timestamp_t t_end_select = get_timestamp();
    double secs_selection = (t_end_select - t_start_select) / 1000000.0L;
    std::cout << "selection stage time " << secs_selection << std::endl;

    double avg_acc = 0;
    for (int p = 0; p < popul.size(); p++) {
      avg_acc += popul[p].score;
    }
    avg_acc = avg_acc / popul.size();

    if (std::abs(avg_acc - popul[best_idx].score) < 0.1) {
      //equal_value = true;
    }

    if (false) {
      // quiero imprimir el score vs rmsd de cada Individuo que me importa para mostrarlos en una grafica.
      double best_rmsd = 0;
      int nearest = calculate_distances_popul->find_nearest(trial_popul[best_idx], popul);
      std::vector<double> rmsd_native_popul = calculate_distances_popul->run_rmsd_to_native(popul);
      std::vector<double> rmsd_native_trial = calculate_distances_popul->run_rmsd_to_native(trial_popul);
      std::cout << "[BEST_IND_ANALYSIS] { best : ( " << (-1 * SCORE_ERROR_FIXED) + popul[best_idx].score << ", " << rmsd_native_popul[best_idx] << " ) , " << " ";
      std::cout << "trial : (" << (-1 * SCORE_ERROR_FIXED) + trial_popul[best_idx].score << " , " << rmsd_native_trial[best_idx] << ") , " << " ";
      std::cout << "parents : [ (" << (-1 * SCORE_ERROR_FIXED) + popul[parents.x1].score <<"," <<rmsd_native_popul[parents.x1] << " ) , "  << "(" << (-1 * SCORE_ERROR_FIXED) + popul[parents.x2].score <<"," <<rmsd_native_popul[parents.x2] << " ), " << "(" << (-1 * SCORE_ERROR_FIXED) + popul[parents.x3].score <<"," <<rmsd_native_popul[parents.x3] << " ) ] , "  << " ";
      std::cout << "nearest : (" << (-1 * SCORE_ERROR_FIXED) + popul[nearest].score << "," << rmsd_native_popul[nearest]  << " ) } " << std::endl;

    }
    print_generation_information(gen_count, new_best_found);
    ++gen_count;

    timestamp_t t_final = get_timestamp();
    double secs_t = (t_final - t_ini) / 1000000.0L;
    std::cout << "gen stage time " << secs_t << std::endl;

  }
}


bool MPIfastCrowdingMoverDE::select_population(const std::vector<Individual>& trial_popul) {
  trial_sucess_n = 0;
  // Crowding DE: Rene Thomsen Algorithm
  // 1. Each trial can be replace only its nearest individual if the energy is better

  bool new_best_found = false;
  double previous_best = popul[best_idx].score;


  std::cout << "[calculated_nearest] ";
  for (int i = 0; i < trial_popul.size(); i++) {
    int target_ind = 0;
    //target_ind = calculate_distances_popul->find_nearest(trial_popul[i], popul);
    target_ind = nearest_inds[i];
    std::cout << target_ind << " ";
    if ( trial_popul[i].score < popul[target_ind].score ) {
      popul[target_ind] = trial_popul[i];
      new_best_found = true;
      trial_sucess_n++;
    }
  }

  std::cout << std::endl;

  avg_acc = 0;
  best_idx = 0;
  best = popul[0].score;
  //  std::cout << "check_popul ";
  for (int i = 0; i < popul.size(); i++) {
    //std::cout << (-1 * SCORE_ERROR_FIXED) + popul[i].score << " ( "<< popul[i].score << " ) ";
    avg_acc += popul[i].score;
    if (popul[i].score < popul[best_idx].score) {
      best_idx = i;
      best = popul[best_idx].score;
      if (best < previous_best) {
  	new_best_found = true;
      }
    }
  }
  //std::cout << std::endl;
  //std::cout << "best " << (-1 * SCORE_ERROR_FIXED) + best << " at " << best_idx << " " << popul[best_idx].score << std::endl;

  copy_popul_fitness_for_print = popul;

  return new_best_found;
}

void
MPISeedsSantosMoverDE::print_current_subpopuls(std::vector<Individual> popul, std::vector<int> current_seed,  std::vector<std::vector<Individual> > subpopuls) {
  std::cout << "[CURRENT_SEEDS] ";
  for (int i = 0; i < current_seed.size(); i++) {
    std::cout << i << " : " <<    (-1 * SCORE_ERROR_FIXED) +  popul[current_seed[i]].score << " , ";
  }
  std::cout << std::endl;

  for (int i = 0; i < subpopuls.size(); i++) {
    std::vector<double> rmsd_subpopul = calculate_distances_popul->run_rmsd_to_native(subpopuls[i]);
    std::cout << "subpopul " << i << " energy [";
    for (int j = 0; j < subpopuls[i].size(); j++) {
      if (j == subpopuls[i].size() - 1) {
	std::cout << " " <<  (-1 * SCORE_ERROR_FIXED) + subpopuls[i][j].score << " ";
      }else {
	std::cout << " " <<  (-1 * SCORE_ERROR_FIXED) + subpopuls[i][j].score << " ,";
      }
    }
    std::cout << "]" << std::endl;
    std::cout << "subpopul " << i << " rmsd [";
    for (int j = 0; j < subpopuls[i].size(); j++) {
      if (j == subpopuls[i].size() - 1) {
	std::cout << " " << rmsd_subpopul[j] << " ";
      }else {
	std::cout << " " << rmsd_subpopul[j] << " ,";
      }
    }
    std::cout << "]" << std::endl;
  }

}

inline  bool sortbysec(const std::pair<int,double> &a,
		       const std::pair<int,double> &b) {
    return (a.second < b.second);
  }


std::vector<int>
MPISeedsSantosKmeansMoverDE::build_kmeans_seeds(std::vector<Individual> popul, std::vector<int> current_seed, std::vector<int> processed_inds, int m_size, int max_subpopuls, std::map<int, std::vector<int> > population_per_seed ) {
  std::vector<int> kmeans_seeds = current_seed;

  for (int i = 0; i < current_seed.size(); i++) {
    std::vector<int> indexes_for_popul = population_per_seed[current_seed[i]];
    std::vector<Individual> temp_popul;
    std::vector<std::pair<int, double> > record_values;
    for (int j = 0; j < indexes_for_popul.size(); j++) {
      temp_popul.push_back(popul[indexes_for_popul[j]]);
    }

    for (int j = 0; j < temp_popul.size(); j++) {
      double dist = 0;
      if (app_options.get<std::string>("Protocol.distance_strategy") == "euclidean_partial_mario") {
	dist = boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->pose_distance_calculation(popul_pdb[ current_seed[i] ], popul_pdb[ indexes_for_popul[j] ] );
      } else {
	dist = calculate_distances_popul->current_distance_calculation(*popul_pdb[current_seed[i]], *popul_pdb[indexes_for_popul[j]] );
      }

      std::pair<int, double> calculate_values = std::pair<int, double>(indexes_for_popul[j], dist);

      record_values.push_back(calculate_values);
    }

    std::sort(record_values.begin(), record_values.end(), sortbysec );
    int new_index_mean = record_values[record_values.size() / 2].first;
    if ( record_values[record_values.size() / 2].second == 0.0) {
      for (int k = record_values.size() / 2; k < record_values.size(); k++) {
	if (  record_values[k].second != 0.0 ) {
	  new_index_mean = record_values[k].first;
	  break;
	}
      }
    }
    kmeans_seeds[i] = new_index_mean;
  }

  return kmeans_seeds;
}

void
MPISeedsSantosMoverDE::build_population_per_seed(std::vector<Individual> popul, std::vector<int> current_seed, std::vector<int> processed_inds, int m_size, int max_subpopuls, std::map<int, std::vector<int> >& population_per_seed ) {
  int min = 1000;
  int min_index = 0;

  for (int i = 0; i < popul.size(); i++) {
    if (processed_inds[i] == 0) {
      min = 1000;
      for (int k = 0; k < current_seed.size(); k++) {
	int subpopul_size = population_per_seed[current_seed[k]].size();
	double dist = 0;
	if (subpopul_size < m_size) {
	  if (app_options.get<std::string>("Protocol.distance_strategy") == "euclidean_partial_mario") {
	    dist = boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->pose_distance_calculation(popul_pdb[current_seed[k]], popul_pdb[i] );
	  } else {
	    dist = calculate_distances_popul->current_distance_calculation(*popul_pdb[current_seed[k]], *popul_pdb[i] );
	  }
	  if (dist < min) {
	    min = dist;
	    min_index = k;
	  }
	}
      }
      population_per_seed[current_seed[min_index]].push_back(i);
    }
  }

  current_seed.resize(max_subpopuls);
}


std::vector<std::vector<Individual> >
MPISeedsSantosMoverDE::create_seeds(std::vector<Individual> popul, std::vector<core::pose::PoseOP> popul_pdb) {
  std::vector<IndIndex> sort_popul_index;
  for (int i = 0; i < popul.size(); i++) {
    sort_popul_index.push_back(IndIndex(i, popul[i]));
  }

  std::vector<Individual> sort_popul = popul;
  std::sort(sort_popul.begin(), sort_popul.end(), popul_comparer );
  PopulSeedsComparer popul_seeds_comparer;
  std::sort(sort_popul_index.begin(), sort_popul_index.end(), popul_seeds_comparer );

  calculate_distances_popul->build_pdb_population(sort_popul, popul_pdb);
  if (app_options.get<std::string>("Protocol.distance_strategy") == "euclidean_partial_mario") {
  boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->build_inter_distances_of_population(sort_popul);
  }

  std::vector<int> current_seed;
  std::vector<int> index_seeds_at_popul;
  current_seed.push_back(sort_popul_index[0].index_);
  index_seeds_at_popul.push_back(sort_popul_index[0].index_);
  //current_cde_radius = 0.03;
  current_cde_radius = fit_rad;
  //  current_cde_radius = 0.01;
  std::map<int, std::vector<int> > population_per_seed;
  std::vector<int> processed_inds(sort_popul.size(), 0);
  processed_inds[sort_popul_index[0].index_] = 1;
  population_per_seed[current_seed[0]] = std::vector<int>();
  population_per_seed[current_seed[0]].push_back( current_seed[0] );
  int max_subpopuls = 8;
  int m_size = popul.size() / max_subpopuls;
  int ind_size = popul[0].vars.size();

  for (int i = 1; i < sort_popul.size(); i++) {
    if (processed_inds[sort_popul_index[i].index_] == 0) {
      bool found = false;
      for (int j = 0; j < current_seed.size() && !found; j++) {
	double dist = 0;
	if (population_per_seed[current_seed[j]].size() < m_size ) {
	  if (app_options.get<std::string>("Protocol.distance_strategy") == "euclidean_partial_mario") {
	    dist = boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->pose_distance_calculation(popul_pdb[current_seed[j]],popul_pdb[i]);
	  } else {
	    dist = calculate_distances_popul->current_distance_calculation(*popul_pdb[current_seed[j]],*popul_pdb[i]);
	  }
	  if ( dist <= current_cde_radius) {
	    found = true;
	    population_per_seed[current_seed[j]].push_back(sort_popul_index[i].index_);
	    processed_inds[sort_popul_index[i].index_]  = 1;
	    break;
	  }
	}
      }

      if (!found) {
	current_seed.push_back(sort_popul_index[i].index_);
	processed_inds[sort_popul_index[i].index_]  = 1;
	index_seeds_at_popul.push_back(sort_popul_index[i].index_);
	population_per_seed[sort_popul_index[i].index_] = std::vector<int>();
	population_per_seed[sort_popul_index[i].index_].push_back(sort_popul_index[i].index_);
      }
    }
  }

  // at this point, only seeds should be marked as processed inds
  // also, index seeds at popul should contain only the index (at popul) of the processed inds (only the current seeds)
  // Limit current seeds to max_subpopuls
  for (int i = 0; i < index_seeds_at_popul.size() && i < max_subpopuls; i++) {
    processed_inds[index_seeds_at_popul[i]] = 0;
  }
  index_seeds_at_popul.resize(max_subpopuls);
  current_seed.resize(max_subpopuls);
  for (int i = 0; i < popul.size(); i++) {
    if (std::find(std::begin(index_seeds_at_popul), std::end(index_seeds_at_popul), i ) != std::end(index_seeds_at_popul) ) {
      processed_inds[i] = 1;
    } else {
      processed_inds[i] = 0;
    }
  }

  for (int i = 0; i < max_subpopuls; i++) {
    population_per_seed[current_seed[i]].resize(1);
  }
  index_seeds_at_popul.resize(max_subpopuls);

  // insert each popul ind at each corresponded popul
  build_population_per_seed(popul, current_seed, processed_inds, m_size, max_subpopuls, population_per_seed );
  std::cout << "print subpopuls" << std::endl;
  std::vector<std::vector<Individual> > temp_subpopuls;
  int count = 0;
  //create temp_subpopuls
  for (int i = 0; i < current_seed.size(); i++) {
    temp_subpopuls.push_back(std::vector<Individual>());
    for (int j = 0; j < population_per_seed[current_seed[i]].size(); j++) {
      temp_subpopuls[i].push_back(popul[population_per_seed[current_seed[i]][j]]);
      count++;
    }
  }
  print_current_subpopuls(popul, current_seed, temp_subpopuls);

  std::vector<std::vector<Individual> > subpopuls = temp_subpopuls;

  //disturb equal inds to seed subpopul
  // for (int i = 0; i < subpopuls.size(); i++) {
  //   for (int j = 1; j < subpopuls[i].size(); j++) {
  //     //if (std::abs(popul[current_seed[i]].score - subpopuls[i][j].score) > 0.01) {
  // 	if (std::abs(subpopuls[i][0].score - subpopuls[i][j].score) > 0.01) {
  // 	init_popul_->disturb_individual(subpopuls[i][j], ind_size);
  //     }
  //   }
  // }

  // score & print subpopuls
  for (int i = 0; i < subpopuls.size(); i++) {
    subpopuls[i] = mpi_calculator->run( subpopuls[i], 1 );
  }

  return subpopuls;
}



std::vector<std::vector<Individual> >
MPISeedsSantosKmeansMoverDE::create_seeds(std::vector<Individual> popul, std::vector<core::pose::PoseOP> popul_pdb) {
  std::vector<IndIndex> sort_popul_index;
  for (int i = 0; i < popul.size(); i++) {
    sort_popul_index.push_back(IndIndex(i, popul[i]));
  }

  std::vector<Individual> sort_popul = popul;
  std::sort(sort_popul.begin(), sort_popul.end(), popul_comparer );
  PopulSeedsComparer popul_seeds_comparer;
  std::sort(sort_popul_index.begin(), sort_popul_index.end(), popul_seeds_comparer );

  calculate_distances_popul->build_pdb_population(sort_popul, popul_pdb);
  if (app_options.get<std::string>("Protocol.distance_strategy") == "euclidean_partial_mario") {
  boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->build_inter_distances_of_population(sort_popul);
  }

  std::vector<int> current_seed;
  std::vector<int> index_seeds_at_popul;
  current_seed.push_back(sort_popul_index[0].index_);
  index_seeds_at_popul.push_back(sort_popul_index[0].index_);
  //current_cde_radius = 0.03;
  current_cde_radius = fit_rad;
  //  current_cde_radius = 0.01;
  std::map<int, std::vector<int> > population_per_seed;
  std::vector<int> processed_inds(sort_popul.size(), 0);
  processed_inds[sort_popul_index[0].index_] = 1;
  population_per_seed[current_seed[0]] = std::vector<int>();
  population_per_seed[current_seed[0]].push_back( current_seed[0] );
  int max_subpopuls = 8;
  int m_size = popul.size() / max_subpopuls;
  int ind_size = popul[0].vars.size();

  for (int i = 1; i < sort_popul.size(); i++) {
    if (processed_inds[sort_popul_index[i].index_] == 0) {
      bool found = false;
      for (int j = 0; j < current_seed.size() && !found; j++) {
	double dist = 0;
	if (population_per_seed[current_seed[j]].size() < m_size ) {
	  if (app_options.get<std::string>("Protocol.distance_strategy") == "euclidean_partial_mario") {
	    dist = boost::dynamic_pointer_cast<CalculateEuclideanMarioPartialDistancePopulation>(calculate_distances_popul)->pose_distance_calculation(popul_pdb[current_seed[j]],popul_pdb[i]);
	  } else {
	    dist = calculate_distances_popul->current_distance_calculation(*popul_pdb[current_seed[j]],*popul_pdb[i]);
	  }
	  if ( dist <= current_cde_radius) {
	    found = true;
	    population_per_seed[current_seed[j]].push_back(sort_popul_index[i].index_);
	    processed_inds[sort_popul_index[i].index_]  = 1;
	    break;
	  }
	}
      }

      if (!found) {
	current_seed.push_back(sort_popul_index[i].index_);
	processed_inds[sort_popul_index[i].index_]  = 1;
	index_seeds_at_popul.push_back(sort_popul_index[i].index_);
	population_per_seed[sort_popul_index[i].index_] = std::vector<int>();
	population_per_seed[sort_popul_index[i].index_].push_back(sort_popul_index[i].index_);
      }
    }
  }

  // at this point, only seeds should be marked as processed inds
  // also, index seeds at popul should contain only the index (at popul) of the processed inds (only the current seeds)
  // Limit current seeds to max_subpopuls
  for (int i = 0; i < index_seeds_at_popul.size() && i < max_subpopuls; i++) {
    processed_inds[index_seeds_at_popul[i]] = 0;
  }
  index_seeds_at_popul.resize(max_subpopuls);
  current_seed.resize(max_subpopuls);
  for (int i = 0; i < popul.size(); i++) {
    if (std::find(std::begin(index_seeds_at_popul), std::end(index_seeds_at_popul), i ) != std::end(index_seeds_at_popul) ) {
      processed_inds[i] = 1;
    } else {
      processed_inds[i] = 0;
    }
  }

  for (int i = 0; i < max_subpopuls; i++) {
    population_per_seed[current_seed[i]].resize(1);
  }
  index_seeds_at_popul.resize(max_subpopuls);

  // insert each popul ind at each corresponded popul
  build_population_per_seed(popul, current_seed, processed_inds, m_size, max_subpopuls, population_per_seed );
  std::cout << "print before kmeans" << std::endl;
  std::vector<std::vector<Individual> > temp_subpopuls;
  int count = 0;
  //create temp_subpopuls
  for (int i = 0; i < current_seed.size(); i++) {
    temp_subpopuls.push_back(std::vector<Individual>());
    for (int j = 0; j < population_per_seed[current_seed[i]].size(); j++) {
      temp_subpopuls[i].push_back(popul[population_per_seed[current_seed[i]][j]]);
      count++;
    }
  }
  print_current_subpopuls(popul, current_seed, temp_subpopuls);

  bool kmeans_strategy = false;
  std::vector<std::vector<Individual> > subpopuls;
  if (kmeans_strategy == true) {
    std::vector<int> seed_kmeans = current_seed;
    std::vector<int> kmeans_processed_inds = processed_inds;
    std::map<int, std::vector<int> > kmeans_popul_per_seed = population_per_seed;
    for(int i = 0; i < 10; i++) {
      seed_kmeans = build_kmeans_seeds(popul, seed_kmeans, kmeans_processed_inds, m_size, max_subpopuls, kmeans_popul_per_seed );

      for (int i = 0; i < kmeans_processed_inds.size(); i++) {
	if (std::find(seed_kmeans.begin(), seed_kmeans.end(), i) == seed_kmeans.end())  {
	  kmeans_processed_inds[i] = 0;
	} else {
	  kmeans_processed_inds[i] = 1;
	}
      }
      kmeans_popul_per_seed = std::map<int, std::vector<int> >()  ;
      for (int i = 0; i < seed_kmeans.size(); i++) {
	kmeans_popul_per_seed[seed_kmeans[i] ] = std::vector<int>() ;
	kmeans_popul_per_seed[seed_kmeans[i] ].push_back(seed_kmeans[i]);
      }

      build_population_per_seed(popul, seed_kmeans, kmeans_processed_inds, m_size, max_subpopuls, kmeans_popul_per_seed);
    }

    population_per_seed = kmeans_popul_per_seed;
    current_seed = seed_kmeans;
    count = 0;
    //create subpopuls
    for (int i = 0; i < current_seed.size(); i++) {
      subpopuls.push_back(std::vector<Individual>());
      for (int j = 0; j < population_per_seed[current_seed[i]].size(); j++) {
	subpopuls[i].push_back(popul[population_per_seed[current_seed[i]][j]]);
	count++;
      }
    }

    for (int i = 0; i < subpopuls.size(); i++) {
      std::vector<Individual> sort_popul = subpopuls[i];
      std::sort(sort_popul.begin(), sort_popul.end(), popul_comparer );
      subpopuls[i] = sort_popul;
    }

    std::cout << "after kmeans" << std::endl;
    print_current_subpopuls(popul, current_seed, subpopuls);
  } else {
    subpopuls = temp_subpopuls;
  }


  //disturb equal inds to seed subpopul
  // for (int i = 0; i < subpopuls.size(); i++) {
  //   for (int j = 1; j < subpopuls[i].size(); j++) {
  //     //if (std::abs(popul[current_seed[i]].score - subpopuls[i][j].score) > 0.01) {
  // 	if (std::abs(subpopuls[i][0].score - subpopuls[i][j].score) > 0.01) {
  // 	init_popul_->disturb_individual(subpopuls[i][j], ind_size);
  //     }
  //   }
  // }

  // score & print subpopuls
  for (int i = 0; i < subpopuls.size(); i++) {
    subpopuls[i] = mpi_calculator->run( subpopuls[i], 1 );
  }

  return subpopuls;
}


#endif
