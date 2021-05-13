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
#include "./helper_apps/FileToPDBPrinter.hpp"

void init_rosetta() {
  std::vector<std::string> arguments = {"./bin/app", "@flags"};
  std::vector<char*> aux_argv;
  for (const auto& arg : arguments) {
   aux_argv.push_back((char*)arg.data());
  }
  aux_argv.push_back(nullptr);

  devel::init(aux_argv.size() - 1, aux_argv.data());
}



void codified_angles_to_pdb(int argc, char** argv) {
  init_rosetta();
  FileToPDBPrinter app_printer;
  app_printer.read(argv[1]);
}


int main(int argc, char *argv[]) {
  std::cout << "init " << std::endl;
  srand ( time(NULL) );

  codified_angles_to_pdb(argc, argv);
  return 0;
}
