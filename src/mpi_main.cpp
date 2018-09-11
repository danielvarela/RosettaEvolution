#include <iostream>
#include <devel/init.hh>
#include <core/pose/Pose.hh>


#include "mpi.h"
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <vector>

#define NUMBER_OF_JOBS 12

namespace mpi = boost::mpi;


#include "PoseFiles/DEoperator.hpp"
#include "PoseFiles/moves/DifferentialEvolutionMover.hpp"
#include "PoseFiles/moves/CalculateRmsdDistancePopulation.hpp"

class MPIDifferentialEvolution {
public:
  mpi::communicator world;

  int slaves_;
  MPIDifferentialEvolution(int slaves) {
    slaves_ = slaves;
  }
  void run(std::vector<int> popul_) {
    // Initialize requests
    unsigned int job_id = 0;
    mpi::request reqs[world.size() - 1];

    std::vector<int> vector_result(popul_.size() * 2);
    // Send initial jobs

    int calc_size = 3;
    std::vector<int> send_popul = popul_;

    int idx_count = 0;
    for (unsigned int dst_rank = 1; dst_rank < world.size(); ++dst_rank) {
      std::cout << "[MASTER] Sending job " << job_id
		<< " to SLAVE " <<  dst_rank << "\n";
      // Send job to dst_rank [nonblocking]
      std::vector<int> values;
      if (idx_count + calc_size == world.size()) {
	values.insert(values.begin(), popul_.begin() + idx_count, popul_.end());
      } else {
	values.insert(values.begin(), popul_.begin() + idx_count, popul_.begin() + idx_count + calc_size);
      }

      std::cout << "values size " << values.size() << std::endl;

      world.isend(dst_rank, 0, dst_rank);

      world.isend(dst_rank, 0, &values.front(), values.size() );
      // Post receive request for new jobs requests by slave [nonblocking]
      std::vector<int>::iterator it_values = vector_result.begin() + idx_count;
      reqs[dst_rank - 1] = world.irecv(dst_rank, dst_rank, *it_values);
      ++job_id;
      idx_count = idx_count + calc_size;
    }

    // Listen for the remaining jobs, and send stop messages on completion.
    // bool all_done = false;
    // while (!all_done) {
    //   all_done = true;
    //   for (unsigned int dst_rank = 1; dst_rank < world.size(); ++dst_rank) {
    // 	if (reqs[dst_rank].test()) {
    // 	  // Tell the slave that it can exit.
    // 	  bool stop = true;
    // 	  world.isend(dst_rank, 0, stop);
    // 	} else {
    // 	  all_done = false;
    // 	}
    //   }
    //   usleep(1000);
    // }
    mpi::wait_all(reqs, reqs + 2);


    for (unsigned int dst_rank = 1; dst_rank < world.size(); ++dst_rank) {
      bool stop = false;
      world.isend(dst_rank, 0, stop);
    }

    for (int i = 0; i < vector_result.size(); i++) {
      std::cout << vector_result[i] << " ";
    }
    std::cout << std::endl;



    std::cout << "[MASTER] Handled all jobs, killed every process.\n";
  }
};

class MPICalculator {
public:

  mpi::communicator world;
  int id_;

  MPICalculator(int rank) {
    id_ = rank;
  }

  void run() {
    //wait for work, calculate, send back the results to master and continue wainting
    std::cout << "thread with rank " << id_ << " is slave " << std::endl;
    bool stop = false;
    while(!stop) {
      // Wait for new jobs
      int tag;
      world.recv(0, 0, tag);

      std::vector<int> vector_in(3);
      world.recv(0, 0, &vector_in.front(), 3);
      std::cout << "[SLAVE: " << world.rank()
		<< "] Received job " << tag << " from MASTER.\n";
      // Perform "job"
      for (int i = 0; i < vector_in.size(); i++) {
	std::cout << vector_in[i] << std::endl;
	vector_in[i] = vector_in[i] + 1;
      }

      // Notify master that the job is done
      world.send(0, tag, &vector_in.front(), vector_in.size());
      // Check if a new job is coming
      world.recv(0, 0, stop);
    }

  }
};


int main(int argc, char** argv) {
  srand ( time(NULL) );
  mpi::environment env(argc, argv);
  mpi::communicator world;

  std::vector<int> popul_(9, 1);
  if (world.rank() == 0) {
    for (int i = 0; i < popul_.size(); i++) {
      popul_[i] = i;
    }
  }

  // int indsPerProcess = popul_.size() / world.size();
  // std::vector<std::vector<int> > populByProcess(world.size(), std::vector<int>());
  // for (int i = 0; i < world.size(); i++) {
  //   for (int k = 0, populIndex = i * indsPerProcess + k; k < indsPerProcess && populIndex < popul_.size() ; k++, populIndex++) {
  //     populByProcess[i].push_back(popul_[populIndex]);
  //   }
  // }

  // std::vector<std::vector<int> > result;
  // scatter(world, populByProcess, result, 0);

  std::map<int, int> map_tags;
  if (world.rank() == 0) {
    std::vector<boost::mpi::request> reqs;
    std::vector<double> num_popul(world.size() * 3, 1);

    int index = 0;
    std::map<int, std::vector<double> > map_populs;
    int array_size;
    for (int i = 1; i < world.size(); i++) {
      std::vector<double> curr_popul = std::vector<double>(num_popul.begin() , num_popul.begin() + 3 );
      map_populs[i] = curr_popul;
      array_size = curr_popul.size();
    }

    for (int i = 1; i < world.size(); i++) {
      std::cout << "send to " << i << " index " << index << " an array of size " << array_size << " " << std::endl;
      map_tags[i] = index;
      // world.send(i, 0, array_size);
      world.isend(i, 0, index);
      world.isend(i, 0, &map_populs[i][0] ,   map_populs[i].size() );
      index = index + (world.size() - 1);
    }

    std::vector<std::vector<double> > res_popul(3, std::vector<double>(3));
    for (int i = 1; i < world.size(); i++) {
      std::cout << "wait for " << i << " index " << map_tags[i] << " array_size " << res_popul[i - 1].size() <<  std::endl;
      //      res_popul[i - 1] = std::vector<double>(3);
      reqs.push_back(world.irecv(i, 0, &res_popul[i - 1][0] , res_popul[i - 1].size()  ) ) ;
    }

    //    mpi::wait_all( std::begin(reqs) , std::end(reqs) );
    bool all_done = false;
    while (!all_done) {
      all_done = true;
      for (unsigned int dst_rank = 1; dst_rank < world.size(); ++dst_rank) {
	//    	  std::cout << "slave " << dst_rank - 1 << " testing " << std::endl;
    	if ( reqs[dst_rank - 1].test() ) {

    	  // Tell the slave that it can exit.
    	  bool stop = false;
    	  //    	  world.isend(dst_rank, 0, stop);
    	} else {
    	  all_done = false;
    	}
      }
      usleep(1000);
    }


    for (int i = 0; i < res_popul.size(); i++) {
      for (std::vector<double>::iterator it = res_popul[i].begin(); it != res_popul[i].end(); ++it) {
	std::cout << *it << "!" << std::endl;
      }

    }

  } else {
    bool stop = true;
    while (stop) {
      int size;
      // world.recv(0, 0, size);
      int index;
      world.recv(0, 0, index);
      std::vector<double> popul(size);
      world.recv(0, 0, &popul[0], popul.size() );

      for (int i = 0; i < popul.size(); i++) {
	popul[i] = world.rank();
      }

      std::cout << "slave " << world.rank() << " index " << index << " sends " << popul.size() << std::endl;
      world.send(0, 0, &popul[0], popul.size() ) ;
      //  world.irecv(0, 0, stop);
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
