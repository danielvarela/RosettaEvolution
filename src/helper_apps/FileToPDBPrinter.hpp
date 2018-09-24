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

   map<string, string> read_configuration_strategy(tokenizer& tokens) {
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

   std::vector<double> read_individual_strategy(tokenizer& tokens) {

    vector<double> vars;
    tokenizer::iterator i = tokens.begin();
    std::advance(i , 1);
    for (  ; i != tokens.end(); i++) {
      vars.push_back(atof((*i).c_str()));
    }

    return vars;
  }

  vector<Individual> population;
  std::map<string, string> configuration;
  boost::shared_ptr<DE_Operator> de;
  boost::shared_ptr<FitFunction> ffxn;
  boost::shared_ptr<PoseFunction> pfunc;

   void read(std::string path) {
    ifstream inFile;
    inFile.open(path.c_str());
    if (!inFile) {
      cerr << "Unable to open file datafile.txt";
      exit(1);   // call system to stop}}
    }
    vector<double> vars;
    std::string line;
    boost::char_separator<char> sep(" ");
    Individual ind;

    while (std::getline(inFile, line)) {
      tokenizer tokens(line, sep);
      if (tokens.begin() != tokens.end()) {
	if (*tokens.begin() == "[IND]") {
          de->init_popul->disturb_individual(ind, de->ffxn->D() );
	  ind.vars = read_individual_strategy(tokens);
	  population.push_back(ind);
	}
	if (*tokens.begin() == "[CONF]") {
	  configuration = read_configuration_strategy(tokens);
	  init_differential_evolution();
	}
      }
    }

    // print_at_screen();
    std::cout << "init de operator " << std::endl;
    std::vector<core::pose::PoseOP> popul_pdb;
    build_pdb_population(popul_pdb);
    // de, population y popul_pdb ->
    // deberia calcular distancias en rmsd, partial_mario y diff absoluto
    // construir json con todos los datos e imprimirlo en un fichero que lee el servidor
    CalculateDistancePopulationPtr distance_calculator = de->use_distances_strategy("rmsd_native_diff");
    //distance_calculator->fit_radius = configuration["fit_radius"];
    distance_calculator->fit_radius = 10.0;
    CalculateDistancePopulation::DistancesResult result = distance_calculator->run(population);


    std::ofstream ofs ( string("/home/dvarela/Code/RosettaEvolution/neighs_example.json") , std::ofstream::out);
    ofs << "{ \"neighs\" : [" << std::endl;
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
    ofs << "}" << std::endl;

    // shared_fitness = result.shared_fitness;
    // rmsd_to_native = result.rmsd_to_native;
    // neighs_per_ind = result.neigh_per_ind;
    // ind_inter_distance = result.distances_of_population;

    // write_population_pdb(de, population);
    std::cout << "finish de operator " << std::endl;
  }

  void build_pdb_population(std::vector<core::pose::PoseOP>& popul_pdb) {
    popul_pdb.resize(0);
    core::pose::PoseOP local_pose;
    core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score3").c_str());
    for (int i = 0; i < population.size(); ++i) {
      local_pose = de->native_pose_->clone();
      pfunc->fill_pose(local_pose, population[i], de->ss);
      popul_pdb.push_back(local_pose);
      std::cout << "score " << i << " : " << (*scorefxn)(*popul_pdb[i]) << std::endl;
    }
    for (int i = 0; i < population.size(); ++i) {
     std::cout << "score " << i << " : " << (*scorefxn)(*popul_pdb[i]) << std::endl;
    }
  }

  void print_at_screen() {
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

  void init_differential_evolution() {
    de = boost::shared_ptr<DE_Operator>(new DE_Operator(configuration["prot"]));
    core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score3").c_str());
    de->frag_opt.scorefxn = scorefxn;
    de->frag_opt.stage_name = std::string("stage4");
    de->frag_mover = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::greedy_search, de->frag_opt);
    ffxn = boost::shared_ptr<FitFunction>( new PoseFragmentFunction(de->pose_, scorefxn, de->ss, de->frag_mover));
    pfunc = boost::dynamic_pointer_cast<PoseFunction >(ffxn);
  }

   void write_population_pdb(boost::shared_ptr<DE_Operator> de,
                            std::vector<Individual>& popul) {
     core::pose::PoseOP local_pose = de->native_pose_->clone();

    for (int i = 0; i < popul.size(); ++i) {
     pfunc->fill_pose(local_pose, popul[i], de->ss);     
     std::string path = "/home/dvarela/Code/test_popul/pOOpul_ind_"+ std::to_string(i) +".pdb";
     local_pose->dump_pdb(path.c_str());
    }
  }


    static void population_to_file_json(boost::property_tree::ptree pt, std::string score_name, int gen, const std::vector<Individual>& popul);

  static void population_to_file(boost::property_tree::ptree pt, std::string score_name, int gen, const std::vector<Individual>& popul) {
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



};
