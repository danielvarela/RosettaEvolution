#ifndef MASTERROSETTACALCULATOR_H
#define MASTERROSETTACALCULATOR_H

#if(MPI_ENABLED)


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

namespace mpi = boost::mpi;

#include "../PoseFiles/moves/CalculateRmsdDistancePopulation.hpp"

class MasterRosettaCalculator {
public:
  mpi::communicator world;
  boost::shared_ptr<FitFunction> scfxn_ind_with_frags;
  boost::shared_ptr<PoseScoreFunction> simple_ffxn;
  explicit MasterRosettaCalculator(boost::shared_ptr<FitFunction> scfxn_ind_with_frags_input)  {
    scfxn_ind_with_frags = scfxn_ind_with_frags_input;
    boost::shared_ptr<PoseFragmentFunction> aux_ffxn = boost::dynamic_pointer_cast<PoseFragmentFunction>(scfxn_ind_with_frags);
    simple_ffxn = boost::shared_ptr<PoseScoreFunction>( new PoseScoreFunction(aux_ffxn->pose_, aux_ffxn->scorefxn, aux_ffxn->ss, aux_ffxn->frag_mover));
#if(USE_CRYO_EM)
    simple_ffxn = boost::shared_ptr<PoseScoreFunction>( new PoseDensityFunction(aux_ffxn->pose_, aux_ffxn->scorefxn, aux_ffxn->ss, aux_ffxn->frag_mover));
#else
    simple_ffxn = boost::shared_ptr<PoseScoreFunction>( new PoseScoreFunction(aux_ffxn->pose_, aux_ffxn->scorefxn, aux_ffxn->ss, aux_ffxn->frag_mover));
#endif



  }

  std::vector<std::vector<IndMPI> >
  SplitVector(const std::vector<Individual>& vec, size_t n)
  {
    std::vector<std::vector<IndMPI> > outVec(n);

    size_t length = vec.size() / n;
    size_t remain = vec.size() % n;

    size_t begin = 0;
    size_t end = 0;
    for (size_t i = 0; i < std::min(n, vec.size()); ++i)
      {
	outVec[i] = std::vector<IndMPI>();
	end += (remain > 0) ? (length + !!(remain--)) : length;
	for (int j = begin; j < end; j++) {
	  outVec[i].push_back(IndMPI(j, vec[j]));
	}
	begin = end;
      }
    return outVec;
  }

  std::vector<Individual>
  run(std::vector<Individual> num_popul, int mode = 2) {
    std::vector<boost::mpi::request> reqs;
    // desde aqui es la funcion calcular popul
    std::string score_name = scfxn_ind_with_frags->name();
    std::string score_func;
    if (score_name == "complete") {
      score_func = "complete";
    } else {
      score_func = std::string(score_name.substr(score_name.rfind(" ") + 1));
    }

    int array_size = num_popul.size() / (world.size());
    //int mode = 2; // for frag insertion
    std::vector<IndMPI> master_vector_for_work;
    std::vector<std::vector<IndMPI> > vector_parts = SplitVector(num_popul, world.size());
    master_vector_for_work = vector_parts[0];
    for (int i = 1; i < world.size(); i++) {
      world.send(i, 0, array_size);
      world.send(i, 0, mode);
      world.send(i, 0, score_func);
      world.send(i, 0, &vector_parts[i].front(), vector_parts[i].size() );
    }
    // prepare res_popul to receive from mpi workers
    std::vector<std::vector<IndMPI> > res_popul(world.size() - 1, std::vector<IndMPI>(array_size));
    for (int i = 1; i < world.size(); i++) {
      reqs.push_back(world.irecv(i, 0, &res_popul[i - 1].front() , array_size ) ) ;
    }
    // score individuals of master
    std::vector<Individual> result_population(num_popul.size());
    for (int i = 0; i < master_vector_for_work.size(); i++) {
      if (mode == 1) {
	simple_ffxn->score(master_vector_for_work[i].ind);
      } else {
	scfxn_ind_with_frags->score(master_vector_for_work[i].ind);
      }
      result_population[master_vector_for_work[i].index] = master_vector_for_work[i].ind;
    }


    mpi::wait_all( std::begin(reqs) , std::end(reqs) );
    // build result_population with received individuals
    for (int i = 0; i < res_popul.size(); i++) {
      for (std::vector<IndMPI>::iterator it = res_popul[i].begin(); it != res_popul[i].end(); ++it) {
	result_population[it->index] = it->ind;
      }
    }
    return result_population;
  }

};

#endif

#endif /* MASTERROSETTACALCULATOR_H */
