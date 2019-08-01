#include <string>
#include <vector>
#include <map>
#include "FileToPDBPrinter.hpp"
#include "Utils.hpp"
#include <boost/algorithm/string/replace.hpp>

using namespace std;

CalculateDistancePopulation::DistancesResult
FileToPDBPrinter::calculate_distances_and_print(std::string calculate_option, CalculateDistancePopulationPtr distance_calculator, std::ofstream& ofs) {
    calculate_option = std::string("euclidean_partial_mario");
    distance_calculator = de->use_distances_strategy(calculate_option);
    distance_calculator->fit_radius = radius_per_option[calculate_option];
    CalculateDistancePopulation::DistancesResult result = print_distances_calculations(distance_calculator, calculate_option, ofs);
    ofs << " }, " << std::endl;
    return result;
    }

void
FileToPDBPrinter::print_values(CalculateDistancePopulation::DistancesResult result , PrintValuesOption print_option, std::ofstream& ofs) {

    std::vector<double> values_to_print = score_values;
    switch (print_option) {
    case print_score_values: {
      values_to_print = score_values;
      ofs << " \"score_values\" : [ ";
      break;
    }
    case print_rmsd_to_native_values: {
	values_to_print = result.rmsd_to_native;
	ofs << " \"rmsd_values\" : [ ";
      break;
    }
default:
      values_to_print = score_values;
      ofs << " \"score_values\" : [ ";
      break;
    }

    for (int i = 0; i < values_to_print.size(); i++) {
      if (i == (values_to_print.size() - 1 )) {
	ofs << values_to_print[i];
      } else {
	ofs << values_to_print[i] << " , ";
      }
    }
    ofs << " ] " << std::endl;
    if (print_score_values) {
      ofs << " , " << std::endl;
    }

}


map<string, string>
FileToPDBPrinter::read_configuration_strategy(tokenizer& tokens) {
  map<string, string> conf;
  tokenizer::iterator i = tokens.begin();
  std::advance(i , 1);
  for (  ; i != tokens.end(); i++) {
    tokenizer::iterator j = i;
    std::advance(j, 1);
    conf[*i] =boost::lexical_cast<unsigned int>(*j);
    std::advance(i, 1);
  }
  return conf;
}

std::vector<double>
FileToPDBPrinter::read_individual_strategy(tokenizer& tokens) {

  vector<double> vars;
  tokenizer::iterator i = tokens.begin();
  std::advance(i , 1);
  for (  ; i != tokens.end(); i++) {
    vars.push_back(atof((*i).c_str()));
  }

  return vars;
}

CalculateDistancePopulation::DistancesResult
FileToPDBPrinter::print_distances_calculations(CalculateDistancePopulationPtr distance_calculator, std::string calculate_option, std::ofstream& ofs) {
  CalculateDistancePopulation::DistancesResult result = distance_calculator->run(population);
  //   ofs << "{ \"neighs_"+calculate_option+"\" : [" << std::endl;
  ofs << " \"neighs_"+calculate_option+"\" : [" << std::endl;
  for (int i = 0; i < result.neigh_per_ind.size(); i++) {
    std::vector<NeighStruct> neighs = result.neigh_per_ind[i];
    ofs << " [ ";
    for (int j = 0; j < neighs.size(); j++) {
      if (j == (neighs.size() - 1)) {
	ofs << "{ \"ind\" : " << neighs[j].index << " , \"dist\" : " << neighs[j].distance << " } ] ";
      } else {
	ofs << "{ \"ind\" : " << neighs[j].index << " , \"dist\" : " << neighs[j].distance << " } , ";
      }
    }

    if (i == (result.neigh_per_ind.size() -1 ) ) {
      ofs << " ] " << std::endl;
    } else {
      ofs << " , " << std::endl;
    }

  }
  //    ofs << "}" << std::endl;
  ofs << " " << std::endl;
  return result;
}


void
FileToPDBPrinter::process_file( ifstream& inFile, PrintPopulationOption popul_option  ) {
  Individual ind;
  std::string line;
  boost::char_separator<char> sep(" ");
  if (popul_option == print_current_population) {
    population.resize(0);
  } else {
    trial_population.resize(0);
  }

  while (std::getline(inFile, line)) {
    tokenizer tokens(line, sep);
    if (tokens.begin() != tokens.end()) {
      if (*tokens.begin() == "[IND]") {
	de->init_popul->disturb_individual(ind, de->ffxn->D() );
	ind.vars = read_individual_strategy(tokens);
	if (popul_option == print_current_population) {
	  population.push_back(ind);
	} else {
	  trial_population.push_back(ind);
	}
      }
      if (*tokens.begin() == "[CONF]") {
	configuration = read_configuration_strategy(tokens);
	init_differential_evolution();
      }
    }
  }
}

void
FileToPDBPrinter::print_data_to_json() {
  std::ofstream ofs ( string("/home/dvarela/Code/RosettaEvolution/neighs_example.json") , std::ofstream::out);
  std::string calculate_option;
  CalculateDistancePopulationPtr distance_calculator;
  ofs << "{ \"neighs_data\" : { ";
  calculate_option = std::string("rmsd_native_diff");
  calculate_distances_and_print(calculate_option, distance_calculator, ofs);

  calculate_option = std::string("rmsd");
  calculate_distances_and_print(calculate_option, distance_calculator, ofs);

  calculate_option = std::string("euclidean_partial_mario");
  CalculateDistancePopulation::DistancesResult result = calculate_distances_and_print(calculate_option, distance_calculator, ofs);

  print_values(result, print_score_values, ofs);
  print_values(result, print_rmsd_to_native_values, ofs);
  /*
    si se pudo procesar trial_population -> calculo sus valores y los imprimo
   */
  ofs << " } " << std::endl;

  ofs.close();
}

void
FileToPDBPrinter::build_pdb_population(std::vector<core::pose::PoseOP>& popul_pdb) {
  popul_pdb.resize(0);
  core::pose::PoseOP local_pose;
  core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score2").c_str());
  score_values.resize(0);
  for (int i = 0; i < population.size(); ++i) {
    local_pose = de->native_pose_->clone();
    pfunc->fill_pose(local_pose, population[i], de->ss);
    popul_pdb.push_back(local_pose);
    double score =  (*scorefxn)(*popul_pdb[i]);
    score_values.push_back( score  );
    std::cout << "score " << i << " : " << score << std::endl;
  }

}

void
FileToPDBPrinter::print_at_screen() {
  std::cout << "CONFG ";
  for (map<string, string>::iterator it = configuration.begin(); it != configuration.end(); ++it) {
    std::cout << it->first << " " << it->second << " ";
  }
  std::cout << std::endl;
  std::vector<double> vars;
  for (size_t j = 0; j < population.size(); j++) {
    vars = population[j].vars;
    std::cout << "IND " << j << " ";
    for (double i : vars) {
      cout << i << " , ";
    }
    cout << endl;
  }
}


void
FileToPDBPrinter::init_differential_evolution() {
  de = boost::shared_ptr<DE_Operator>(new DE_Operator(configuration["prot"]));
  //  de = boost::shared_ptr<DE_Operator>(new DE_Operator("1elwA"));
  core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score3").c_str());
  de->frag_opt.scorefxn = scorefxn;
  de->frag_opt.stage_name = std::string("stage4");
  de->frag_mover = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::greedy_search, de->frag_opt);
  ffxn = boost::shared_ptr<FitFunction>( new PoseFragmentFunction(de->pose_, scorefxn, de->ss, de->frag_mover));
  pfunc = boost::dynamic_pointer_cast<PoseFunction >(ffxn);
}


void
FileToPDBPrinter::write_population_pdb(boost::shared_ptr<DE_Operator> de,
				       std::vector<Individual>& popul) {
  core::pose::PoseOP local_pose = de->native_pose_->clone();

  for (int i = 0; i < popul.size(); ++i) {
    pfunc->fill_pose(local_pose, popul[i], de->ss);
    std::string path = "/home/dvarela/Code/test_popul/pOOpul_ind_"+ std::to_string(i) +".pdb";
    local_pose->dump_pdb(path.c_str());
  }
}




void
FileToPDBPrinter::population_to_file_json(boost::property_tree::ptree pt, std::string score_name, int gen, const std::vector<Individual>& popul) {
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

   const int dir_err = mkdir( std::string("/home/dvarela/Code/RosettaEvolution/"+output_folder).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (-1 == dir_err) {
      printf("Error creating directory!n");
      exit(1);
    }
 //  if (dirExists( std::string("./"+output_folder).c_str() ) == 0) {
 // }
  std::ofstream ofs ( string("./"+output_folder+"/population_"+sc_name+"_"+ to_string(gen) +".log") , std::ofstream::out);
  //ofs << "[CONF] " << "NP " << NP << " D " << D << " score " << sc_name << " fit_rad "<<  fit_rad << " prot " << prot << " gen " << gen << std::endl;
  ofs << "{" << std::endl;
  ofs << " conf : { " << "NP : " << NP << " ,  D : " << D << " ,  score : " << sc_name << " ,  fit_rad : "<<  fit_rad << " , prot : " << prot << " ,  gen : " << gen << "} , " << std::endl;

  ofs << " individuals : [";
  for (size_t i = 0; i < popul.size(); i++) {
    ofs << " [ ";
    std::vector<double> vars = popul[i].vars;
    for (size_t j = 0; j < vars.size(); j++) {
      if (j == (vars.size() - 1)) {
	ofs << vars[j] << " ] ";
      } else {
	ofs << vars[j] << " , ";
      }
    }
    if (i == (popul.size() - 1)) {
      ofs << std::endl << " ] " << std::endl;
    } else {
      ofs << " , " << std::endl;
    }
  }
  ofs << "}" << std::endl;
  ofs.close();

}
