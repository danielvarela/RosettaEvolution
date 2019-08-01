// Copyright 2018 by vare
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <iterator>
#include <boost/tokenizer.hpp>
#include <boost/property_tree/ptree.hpp>

#include "../Algorithm/DE_types.hpp"
#include "../Controller/DE_Operator.hpp"

#include <boost/algorithm/string/replace.hpp>
using namespace std;

class FileToPDBPrinter {
public:
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  vector<Individual> population, trial_population;
  std::map<string, string> configuration;
  boost::shared_ptr<DE_Operator> de;
  boost::shared_ptr<FitFunction> ffxn;
  boost::shared_ptr<PoseFunction> pfunc;
  std::vector<double> score_values;
  double partial_mario_rad;
  double rmsd_rad;
  double diff_abs_rad;
  std::map<std::string, double> radius_per_option;

  enum PrintPopulationOption {
    print_current_population, print_trial_population
  };

  enum PrintValuesOption {
    print_score_values, print_rmsd_to_native_values
  };

  FileToPDBPrinter() {
  partial_mario_rad = 1.0;
  rmsd_rad = 1.0;
  diff_abs_rad = 1.0;
  char* tmp = getenv( "diff_abs_rad" );

  if ( tmp == NULL ) {
    diff_abs_rad = 1.0;
  } else {
    diff_abs_rad = std::atof(std::string(tmp).c_str());
  }
  tmp = getenv( "rmsd_rad" );

  if ( tmp == NULL ) {
    rmsd_rad = 1.0;
  } else {
    rmsd_rad = std::atof(std::string(tmp).c_str());
  }
  tmp = getenv( "partial_mario_rad" );

  if ( tmp == NULL ) {
    partial_mario_rad = 1.0;
  } else {
    partial_mario_rad = std::atof(std::string(tmp).c_str());
  }


  radius_per_option["diff_abs_rad"] = diff_abs_rad;
  radius_per_option["rmsd_rad"] = rmsd_rad;
  radius_per_option["partial_mario_rad"] = partial_mario_rad;
}
  map<string, string> read_configuration_strategy(tokenizer& tokens);

  std::vector<double> read_individual_strategy(tokenizer& tokens);

  void
  process_file( ifstream& inFile, PrintPopulationOption popul_option  ) ;

  CalculateDistancePopulation::DistancesResult
  calculate_distances_and_print(std::string calculate_option, CalculateDistancePopulationPtr distance_calculator, std::ofstream& ofs);

  void
  print_values(CalculateDistancePopulation::DistancesResult result , PrintValuesOption print_option , std::ofstream& ofs);

  void print_data_to_json();

  void read(std::string path) {
  ifstream inFile;
  inFile.open(path.c_str());
  if (!inFile) {
    cerr << "Unable to open file datafile.txt";
    exit(1);   // call system to stop}}
  }
  process_file(inFile, print_current_population);
  std::string trial_popul_path = path;
  boost::replace_all(trial_popul_path, "population", "trial_population");
  inFile.open(path.c_str());
  if (inFile) {
    process_file(inFile, print_trial_population);
  }

  // print_at_screen();
  std::cout << "init de operator " << std::endl;
  std::vector<core::pose::PoseOP> popul_pdb;
  build_pdb_population(popul_pdb);
  print_data_to_json();
  write_population_pdb(de, population);
  std::cout << "finish de operator " << std::endl;
  }

  CalculateDistancePopulation::DistancesResult
  print_distances_calculations(CalculateDistancePopulationPtr distance_calculator, std::string calculate_option, std::ofstream& ofs);

  void build_pdb_population(std::vector<core::pose::PoseOP>& popul_pdb);

  void print_at_screen();

  void init_differential_evolution();

  void write_population_pdb(boost::shared_ptr<DE_Operator> de,
			    std::vector<Individual>& popul);

  static void population_to_file_json(boost::property_tree::ptree pt, std::string score_name, int gen, const std::vector<Individual>& popul);

  static void population_to_file(boost::property_tree::ptree pt, std::string score_name, int gen, const std::vector<Individual>& popul) {
    //    population_to_file_json(pt, score_name, gen, popul);
    population_to_file(pt, score_name, gen, popul, print_current_population);
  }

  static void population_to_file(boost::property_tree::ptree pt, std::string score_name, int gen, const std::vector<Individual>& popul, PrintPopulationOption print_option) {
    //    population_to_file_json(pt, score_name, gen, popul);
    int NP = pt.get<int>("DE.NP");
    int D = popul[0].vars.size();
    boost::char_separator<char> sep(" ");
    tokenizer tokens(score_name, sep);
    std::string sc_name;
    tokenizer::iterator it = tokens.begin();
    for (; it != tokens.end(); ++it) {
      sc_name = *it;
    }

    double fit_rad = pt.get<double>("Extra.fitrad");
    std::string prot = pt.get<std::string>("Protocol.prot");
    std::string output_folder = pt.get<std::string>("PrintPopul.output");
    if (dirExists( std::string("./"+output_folder).c_str() ) == 0) {
      const int dir_err = mkdir( std::string("./"+output_folder).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      if (-1 == dir_err) {
	printf("Error creating directory!n");
	exit(1);
      }

    }

    std::string file_name = std::string("population");
    if (print_option == print_trial_population) {
      file_name = std::string("trial_population");
    }

    std::ofstream ofs ( string("./"+output_folder+"/"+file_name+"_"+sc_name+"_"+ to_string(gen) +".log") , std::ofstream::out);
    ofs << "[CONF] " << "NP " << NP << " D " << D << " score " << sc_name << " fit_rad "<<  fit_rad << " prot " << prot << " gen " << gen << std::endl;
    for (size_t i = 0; i < popul.size(); i++) {
      ofs << "[IND] ";
      std::vector<double> vars = popul[i].vars;
      for (size_t j = 0; j < vars.size(); j++) {
    	ofs << vars[j] << " ";
      }
      ofs << std::endl;
    }
    ofs.close();
  }



};
