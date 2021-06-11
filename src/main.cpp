#include <iostream>
#include <devel/init.hh>

#include "mpi.h"
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <vector>
#include <map>
#include <string>
#include <core/pose/Pose.hh>
namespace mpi = boost::mpi;

#include "Controller/DE_Operator.hpp"
#include "MpiFiles/UtilsMPI.hpp"
#include "MpiFiles/WorkerProcess.hpp"

#define BOOST_DATE_TIME_NO_LIB
#define NUMBER_OF_JOBS 12


// Includes
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/local_time_adjustor.hpp>
#include <boost/date_time/c_local_time_adjustor.hpp>
#include <boost/optional/optional.hpp>

void
init_rosetta(int res) {
  std::string init_flags = "@flags";
  switch (res) {
  case 2: {
    init_flags = "@flags_2"; break;
  }
  case 5: {
    init_flags = "@flags_5"; break;
  }
  case 10: {
    init_flags = "@flags_10"; break;
  }
  default:
    init_flags = "@flags"; break;
  }
  std::vector<std::string> arguments = {"./bin/mpi", init_flags};
  std::vector<char*> aux_argv;
  for (const auto& arg : arguments) {
   aux_argv.push_back((char*)arg.data());
  }
  aux_argv.push_back(nullptr);

  devel::init(aux_argv.size() - 1, aux_argv.data());
}

void show_options_file_in_log(std::string file_name) {
  std::string line;
  std::cout << "[OPTIONS_FILE_USED]" << std::endl;
  std::ifstream myfile(file_name);
  if (myfile.is_open()) {
    while ( std::getline (myfile,line) ) {
      if (line.size() > 0) {
	  if (line[0] != '#') std::cout << line << '\n';
        }
      }
      myfile.close();
    }
  std::cout << "[END_OPTIONS_FILE]" << std::endl;
}


boost::shared_ptr<DE_Operator>
run_operator(int argc, char** argv) {
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(std::string(argv[argc - 1]), pt);
  std::string prot_name = pt.get<std::string>("Protocol.prot");
  int map_res = pt.get<int>("Protocol.map_res");
  init_rosetta(map_res);
  return boost::shared_ptr<DE_Operator>(new DE_Operator(prot_name, pt) );
}

int
main(int argc, char** argv) {
  using boost::posix_time::ptime;
  using boost::posix_time::second_clock;
  using boost::posix_time::to_simple_string;
  using boost::gregorian::day_clock;
  srand ( time(NULL) );
  mpi::environment env(argc, argv);
  mpi::communicator world;


  if (world.rank() == 0) {
    boost::shared_ptr<DE_Operator> app_operator = run_operator(argc, argv);
    show_options_file_in_log(argv[argc - 1]);
    app_operator->run();
    for (int i = 1; i < world.size(); i++) {
      world.send(i, 0, STOP_TAG);
    }
    ptime todayUtc(day_clock::universal_day(), second_clock::universal_time().time_of_day());
    std::cout << "time " << to_simple_string(todayUtc) << std::endl;
  } else {
    boost::shared_ptr<DE_Operator> app_operator = run_operator(argc, argv);
    //WorkerProcess(app_operator).run();
    WorkerProcessScatterGather(app_operator).run();
    //WorkerProcessEvaluateAndNearest(app_operator).run();
  }
}
