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

#include "../PoseFiles/moves/DE_types.hpp"
#include "../PoseFiles/DEoperator.hpp"

using namespace std;

class FileToPDBPrinter {
public:
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

  FileToPDBPrinter() {

  }

  static void population_to_file(boost::property_tree::ptree pt, std::string score_name, int gen, const std::vector<Individual>& popul) {
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
    
    std::ofstream ofs ( string("./"+output_folder+"/population_"+sc_name+"_"+ to_string(gen) +".log") , std::ofstream::out);
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
    boost::shared_ptr<DE_Operator> de; 
    while (std::getline(inFile, line)) {
      tokenizer tokens(line, sep);
      if (tokens.begin() != tokens.end()) {
	if (*tokens.begin() == "[IND]") {
	  Individual ind;
          de->init_popul->disturb_individual(ind, de->ffxn->D() );
	  ind.vars = individual_strategy(tokens);
	  population.push_back(ind);
	}

	if (*tokens.begin() == "[CONF]") {
	  configuration = configuration_strategy(tokens);
          de = boost::shared_ptr<DE_Operator>(new DE_Operator(configuration["prot"]));
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

    write_population_pdb(de, population);
    std::cout << "finish de operator " << std::endl;
  }

  static void write_population_pdb(boost::shared_ptr<DE_Operator> de,
                            std::vector<Individual>& popul) {
    core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score3").c_str());
    de->frag_opt.scorefxn = scorefxn;
    de->frag_opt.stage_name = std::string("stage4");
    de->frag_mover = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::greedy_search, de->frag_opt);
    boost::shared_ptr<FitFunction> ffxn = boost::shared_ptr<FitFunction>( new PoseFragmentFunction(de->pose_, scorefxn, de->ss, de->frag_mover));
    boost::shared_ptr<PoseFunction> pfunc = boost::dynamic_pointer_cast<PoseFunction >(ffxn);
    core::pose::PoseOP local_pose = de->native_pose_->clone(); 

    for (int i = 0; i < popul.size(); ++i) {
     pfunc->fill_pose(local_pose, popul[i], de->ss);     
     std::string path = "/home/dani/Workspace/RosettaEvolution/popul_pdb/popul_ind_"+ std::to_string(i) +".pdb";
     local_pose->dump_pdb(path.c_str());
    }
  }
};
