#include <iostream>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <vector>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "PoseFiles/DEoperator.hpp"
#include "PoseFiles/moves/DifferentialEvolutionMover.hpp"
#include "PoseFiles/moves/CalculateRmsdDistancePopulation.hpp"

void init_rosetta() {
  std::vector<std::string> arguments = {"./bin/app", "@flags"};
  std::vector<char*> aux_argv;
  for (const auto& arg : arguments) {
   aux_argv.push_back((char*)arg.data());
  }
  aux_argv.push_back(nullptr);

  devel::init(aux_argv.size() - 1, aux_argv.data());
}

void run_operator(int argc, char** argv) {
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(std::string(argv[1]), pt);

  std::cout << pt.get<std::string>("Protocol.stages") << std::endl;
  std::cout << pt.get<std::string>("Protocol.distance_strategy") << std::endl;
  init_rosetta();
  std::string prot_name = pt.get<std::string>("Protocol.prot");
  DE_Operator my_app(prot_name, pt);

  std::cout << " Test differential Evolution Riken " << std::endl;
  std::cout << "protname " << prot_name << std::endl;
  std::cout << "protocol " << pt.get<std::string>("Protocol.name") << std::endl;
  my_app.run();
}

#include "./helper_apps/FileToPDBPrinter.hpp"

void codified_angles_to_pdb(int argc, char** argv) {
  init_rosetta();
  FileToPDBPrinter::print(argv[1]);
}

void run_score(int argc, char** argv) {
 //  devel::init(argc, argv);
  init_rosetta();
  std::cout << " Score protein " << std::string(argv[2]) << std::endl;
  DE_Operator my_app = DE_Operator(std::string(argv[2]));

  //  my_app.init_files(my_app.prot_selection[std::string(argv[2])]);

  core::pose::PoseOP native_ = my_app.get_native_pose();

  core::pose::Pose input_pdb;
  read_pose(std::string(argv[3]), input_pdb);

  core::scoring::ScoreFunctionOP scorefxn3 = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score3").c_str());

  std::cout << "score3 for input_pdb is " << (*scorefxn3)(input_pdb) << std::endl;
  std::cout << "rmsd for input_pdb is " << core::scoring::CA_rmsd(input_pdb, *native_) << std::endl;

}

int main(int argc, char *argv[]) {
  bool score = false;
  std::cout << "init " << std::endl;
  srand ( time(NULL) );

  //  codified_angles_to_pdb(argc, argv);

  if (score) {
    run_score(argc, argv);
  } else {
    run_operator(argc, argv);
  }

  return 0;
}
