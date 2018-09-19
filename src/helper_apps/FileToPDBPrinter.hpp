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

#include "../PoseFiles/moves/DE_types.hpp"
#include "../PoseFiles/DEoperator.hpp"

using namespace std;

class FileToPDBPrinter {
public:
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

  FileToPDBPrinter() {

  }

  static void population_to_file( int NP, int D, std::string score_name, double fit_rad, const std::string & prot, int gen, const std::vector<Individual>& popul) {
    std::ofstream ofs ( string("./temp_4/population"+score_name+"_"+ to_string(gen) +".log") , std::ofstream::out);
    ofs << "[CONF] " << "NP " << NP << " D " << D << " score " << score_name << " fit_rad "<<  fit_rad << " prot " << prot << " gen " << gen << std::endl;
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

  static map<string, string> configuration_strategy(tokenizer& tokens) {

    map<string, string> conf;
    tokenizer::iterator i = tokens.begin();
    advance(i , 1);
    for (  ; i != tokens.end(); i++) {
      tokenizer::iterator j = i;
      advance(j, 1);
      conf[*i] = *(j);
      advance(i, 1);
    }

    return conf;
  }



  static std::vector<double> individual_strategy(tokenizer& tokens) {

    vector<double> vars;
    tokenizer::iterator i = tokens.begin();
    std::advance(i , 1);
    for (  ; i != tokens.end(); i++) {
      vars.push_back(atof((*i).c_str()));
    }

    return vars;
  }


  static void print(std::string path) {
    ifstream inFile;
    inFile.open(path.c_str());
    if (!inFile) {
      cerr << "Unable to open file datafile.txt";
      exit(1);   // call system to stop}}
    }
    vector<double> vars;
    vector<Individual> population;
    std::string line;
    boost::char_separator<char> sep(" ");

    std::map<string, string> configuration;

    while (std::getline(inFile, line)) {
      tokenizer tokens(line, sep);
      if (tokens.begin() != tokens.end()) {
	if (*tokens.begin() == "[IND]") {
	  Individual ind;
	  ind.vars = individual_strategy(tokens);
	  population.push_back(ind);
	}

	if (*tokens.begin() == "[CONF]") {
	  configuration = configuration_strategy(tokens);
	}
      }
    }

    std::cout << "CONFG ";
    for (map<string, string>::iterator it = configuration.begin(); it != configuration.end(); ++it) {
      std::cout << it->first << " " << it->second << " ";
    }
    std::cout << std::endl;
    for (size_t j = 0; j < population.size(); j++) {
      vars = population[j].vars;
      std::cout << "IND " << j << " ";
      for (double i : vars) {
	cout << i << " , ";
      }
      cout << endl;
    }

    std::cout << "init de operator " << std::endl;
    DE_Operator de(configuration["prot"]);
    std::cout << "finish de operator " << std::endl;

  }
};
