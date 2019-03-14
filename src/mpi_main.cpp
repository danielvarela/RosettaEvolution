#include <iostream>
#include <devel/init.hh>
#include <core/pose/Pose.hh>


#include "mpi.h"
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <vector>
#include <map>
#include <string>

#define NUMBER_OF_JOBS 12

namespace mpi = boost::mpi;


#include "PoseFiles/DEoperator.hpp"
#include "PoseFiles/moves/DifferentialEvolutionMover.hpp"
#include "PoseFiles/moves/CalculateRmsdDistancePopulation.hpp"
#include "./mpi_files/UtilsMPI.hpp"

#define BOOST_DATE_TIME_NO_LIB

// Includes
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/local_time_adjustor.hpp>
#include <boost/date_time/c_local_time_adjustor.hpp>
#include <boost/optional/optional.hpp>


void
init_rosetta() {
  std::vector<std::string> arguments = {"./bin/mpi", "@flags"};
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
  //  show_options_file_in_log(argv[argc - 1]);

  std::cout << pt.get<std::string>("Protocol.stages") << std::endl;
  std::cout << pt.get<std::string>("Protocol.distance_strategy") << std::endl;
  std::string prot_name = pt.get<std::string>("Protocol.prot");
  return boost::shared_ptr<DE_Operator>(new DE_Operator(prot_name, pt) );
}

int
main(int argc, char** argv) {
  srand ( time(NULL) );
  mpi::environment env(argc, argv);
  mpi::communicator world;

  init_rosetta();

  std::map<int, int> map_tags;
  if (world.rank() == 0) {
    int stop_tag = -42;
    boost::shared_ptr<DE_Operator> de = run_operator(argc, argv);
    show_options_file_in_log(argv[argc - 1]);
    de->run();
    for (int i = 1; i < world.size(); i++) {
      world.send(i, 0, stop_tag);
    }
    using boost::posix_time::ptime;
    using boost::posix_time::second_clock;
    using boost::posix_time::to_simple_string;
    using boost::gregorian::day_clock;

    ptime todayUtc(day_clock::universal_day(), second_clock::universal_time().time_of_day());
    std::cout << " time " << to_simple_string(todayUtc) << std::endl;


  } else {
    bool stop = true;

#if(USE_CRYO_EM)
    std::string frag_type = "hybrid_mover";
#else
    std::string frag_type = "stage_rosetta_mover";
#endif
    boost::shared_ptr<DE_Operator> de = run_operator(argc, argv);

    boost::shared_ptr<FitFunction> scfxn_ind_with_frags_1,  scfxn_ind_simple_1;
    boost::shared_ptr<FitFunction> scfxn_ind_with_frags_2,  scfxn_ind_simple_2;
    boost::shared_ptr<FitFunction> scfxn_ind_with_frags_3,  scfxn_ind_simple_3;
    boost::shared_ptr<FitFunction> scfxn_ind_with_frags_complete,  scfxn_ind_simple_complete;
    init_functions("stage2" ,frag_type, de, scfxn_ind_with_frags_1, scfxn_ind_simple_1 );
    init_functions("stage3" ,frag_type, de, scfxn_ind_with_frags_2, scfxn_ind_simple_2 );
    init_functions("stage4" ,frag_type, de, scfxn_ind_with_frags_3, scfxn_ind_simple_3 );
    init_functions("complete" ,frag_type, de, scfxn_ind_with_frags_complete, scfxn_ind_simple_complete );

    while (stop) {
      int size;
      world.recv(0, 0, size);
      if (size == -42) {
	break;
      }
      int mode;
      world.recv(0, 0, mode);
      std::string score_func;
      world.recv(0, 0, score_func);
      std::vector<IndMPI> popul(size);
      world.recv(0, 0, &popul.front(), popul.size() );
      ScoreStrategy score_ind( scfxn_ind_with_frags_1,  scfxn_ind_simple_1,
			       scfxn_ind_with_frags_2,  scfxn_ind_simple_2,
			       scfxn_ind_with_frags_3,  scfxn_ind_simple_3,
			       scfxn_ind_with_frags_complete,  scfxn_ind_simple_complete,
			       mode, score_func);
      for (int i = 0; i < popul.size(); i++) {
	score_ind.apply(popul[i].ind);
      }
      //std::cout << "slave " << world.rank() << " sends " << popul.size() << std::endl;
      world.send(0, 0, &popul.front(), popul.size() ) ;
      //      world.recv(0, 0, stop);
    }

 }

}

// OLD MAIN
  // if (world.rank() == 0) {
  //   std::cout << "MASTER NODE " << std::endl;
  //   std::vector<int> popul_(9, 1);
  //   for (int i = 0; i < popul_.size(); i++) {
  //     popul_[i] = i;
  //   }
  //   MPIDifferentialEvolution de(world.size() - 1);
  //   de.run(popul_);
  // } else {
  //   MPICalculator calc(world.rank());
  //   calc.run();
  // }




// OLD MASTER

    // // Initialize requests
    // unsigned int job_id = 0;
    // mpi::request reqs[world.size() - 1];

    // std::vector<int> vector_result(popul_.size() * 2);
    // // Send initial jobs

    // int calc_size = 3;
    // std::vector<int> send_popul = popul_;

    // int idx_count = 0;
    // for (unsigned int dst_rank = 1; dst_rank < world.size(); ++dst_rank) {
    //   std::cout << "[MASTER] Sending job " << job_id
    // 		<< " to SLAVE " <<  dst_rank << "\n";
    //   // Send job to dst_rank [nonblocking]
    //   std::vector<int> values;
    //   if (idx_count + calc_size == world.size()) {
    // 	values.insert(values.begin(), popul_.begin() + idx_count, popul_.end());
    //   } else {
    // 	values.insert(values.begin(), popul_.begin() + idx_count, popul_.begin() + idx_count + calc_size);
    //   }

    //   std::cout << "values size " << values.size() << std::endl;

    //   world.isend(dst_rank, 0, dst_rank);

    //   world.isend(dst_rank, 0, &values.front(), values.size() );
    //   // Post receive request for new jobs requests by slave [nonblocking]
    //   std::vector<int>::iterator it_values = vector_result.begin() + idx_count;
    //   reqs[dst_rank - 1] = world.irecv(dst_rank, dst_rank, *it_values);
    //   ++job_id;
    //   idx_count = idx_count + calc_size;
    // }

    // // Listen for the remaining jobs, and send stop messages on completion.
    // // bool all_done = false;
    // // while (!all_done) {
    // //   all_done = true;
    // //   for (unsigned int dst_rank = 1; dst_rank < world.size(); ++dst_rank) {
    // // 	if (reqs[dst_rank].test()) {
    // // 	  // Tell the slave that it can exit.
    // // 	  bool stop = true;
    // // 	  world.isend(dst_rank, 0, stop);
    // // 	} else {
    // // 	  all_done = false;
    // // 	}
    // //   }
    // //   usleep(1000);
    // // }
    // mpi::wait_all(reqs, reqs + 2);


    // for (unsigned int dst_rank = 1; dst_rank < world.size(); ++dst_rank) {
    //   bool stop = false;
    //   world.isend(dst_rank, 0, stop);
    // }

    // for (int i = 0; i < vector_result.size(); i++) {
    //   std::cout << vector_result[i] << " ";
    // }
    // std::cout << std::endl;



    // std::cout << "[MASTER] Handled all jobs, killed every process.\n";


// GATHER EXAMPLE NOT WORKING RIGHT

  // int start = 0;
  // int slave_size = popul_.size() / (world.size() - 1 );
  // int end = start + slave_size;
  // if (world.rank() == 0) {
  //   std::cout << "MASTER NODE " << std::endl;
  //   std::vector<int> output_numbers;
  //   gather(world, &popul_.front(), popul_.size(), output_numbers, 0);
  //   for (int i = 0; i < output_numbers.size(); i++) {
  //     std::cout << output_numbers[i] << std::endl;
  //   }
  // } else {
  //   start = (world.rank() - 1) * slave_size;
  //   end = start + slave_size;
  //   std::vector<double> new_popul;
  //   new_popul.insert(new_popul.end(), popul_.begin() + start, popul_.begin() + end);
  //   for (int i = 0; i < new_popul.size(); i++) {
  //     new_popul[i] = new_popul[i] * 10;
  //     //      std::cout << new_popul[i] << std::endl;
  //   }
  //   gather(world, &new_popul.front(), slave_size, 0);
  // }

// class MPIDifferentialEvolution {
// public:
//   mpi::communicator world;

//   int slaves_;
//   MPIDifferentialEvolution(int slaves) {
//     slaves_ = slaves;
//   }
//   void run(std::vector<int> popul_) {
//     // Initialize requests
//     unsigned int job_id = 0;
//     mpi::request reqs[world.size() - 1];

//     std::vector<int> vector_result(popul_.size() * 2);
//     // Send initial jobs

//     int calc_size = 3;
//     std::vector<int> send_popul = popul_;

//     int idx_count = 0;
//     for (unsigned int dst_rank = 1; dst_rank < world.size(); ++dst_rank) {
//       std::cout << "[MASTER] Sending job " << job_id
// 		<< " to SLAVE " <<  dst_rank << "\n";
//       // Send job to dst_rank [nonblocking]
//       std::vector<int> values;
//       if (idx_count + calc_size == world.size()) {
// 	values.insert(values.begin(), popul_.begin() + idx_count, popul_.end());
//       } else {
// 	values.insert(values.begin(), popul_.begin() + idx_count, popul_.begin() + idx_count + calc_size);
//       }

//       std::cout << "values size " << values.size() << std::endl;

//       world.isend(dst_rank, 0, dst_rank);

//       world.isend(dst_rank, 0, &values.front(), values.size() );
//       // Post receive request for new jobs requests by slave [nonblocking]
//       std::vector<int>::iterator it_values = vector_result.begin() + idx_count;
//       reqs[dst_rank - 1] = world.irecv(dst_rank, dst_rank, *it_values);
//       ++job_id;
//       idx_count = idx_count + calc_size;
//     }

//     // Listen for the remaining jobs, and send stop messages on completion.
//     // bool all_done = false;
//     // while (!all_done) {
//     //   all_done = true;
//     //   for (unsigned int dst_rank = 1; dst_rank < world.size(); ++dst_rank) {
//     // 	if (reqs[dst_rank].test()) {
//     // 	  // Tell the slave that it can exit.
//     // 	  bool stop = true;
//     // 	  world.isend(dst_rank, 0, stop);
//     // 	} else {
//     // 	  all_done = false;
//     // 	}
//     //   }
//     //   usleep(1000);
//     // }
//     mpi::wait_all(reqs, reqs + 2);


//     for (unsigned int dst_rank = 1; dst_rank < world.size(); ++dst_rank) {
//       bool stop = false;
//       world.isend(dst_rank, 0, stop);
//     }

//     for (int i = 0; i < vector_result.size(); i++) {
//       std::cout << vector_result[i] << " ";
//     }
//     std::cout << std::endl;



//     std::cout << "[MASTER] Handled all jobs, killed every process.\n";
//   }
// };

// class MPICalculator {
// public:

//   mpi::communicator world;
//   int id_;

//   explicit MPICalculator(int rank) {
//     id_ = rank;
//   }

//   void run() {
//     //wait for work, calculate, send back the results to master and continue wainting
//     std::cout << "thread with rank " << id_ << " is slave " << std::endl;
//     bool stop = false;
//     while(!stop) {
//       // Wait for new jobs
//       int tag;
//       world.recv(0, 0, tag);

//       std::vector<int> vector_in(3);
//       world.recv(0, 0, &vector_in.front(), 3);
//       std::cout << "[SLAVE: " << world.rank()
// 		<< "] Received job " << tag << " from MASTER.\n";
//       // Perform "job"
//       for (int i = 0; i < vector_in.size(); i++) {
// 	std::cout << vector_in[i] << std::endl;
// 	vector_in[i] = vector_in[i] + 1;
//       }
//       // Notify master that the job is done
//       world.send(0, tag, &vector_in.front(), vector_in.size());
//       // Check if a new job is coming
//       world.recv(0, 0, stop);
//     }
//   }
// };
