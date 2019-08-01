
#ifndef UTILSMPI_H
#define UTILSMPI_H

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
#define STOP_TAG -42
#define NUMBER_OF_JOBS 12

namespace mpi = boost::mpi;


#include "../Controller/DE_Operator.hpp"
#include "../Controller/FitnessFunctionInitializer.hpp"
#include "../Algorithm/DifferentialEvolutionMover.hpp"
#include "../Movers/CalculateRmsdDistancePopulation.hpp"

inline void init_functions(std::string stage_name, std::string frag_type, boost::shared_ptr<DE_Operator> de, boost::shared_ptr<FitFunction>& scfxn_ind_with_frags, boost::shared_ptr<FitFunction>& scfxn_ind_simple) {
  std::string score_fxn_name = "score3";
  if (stage_name == "stage2") {
    score_fxn_name = "score1";
  }
  if (stage_name == "stage3") {
    score_fxn_name = "score2";
  }
  if (stage_name == "stage4") {
    score_fxn_name = "score3";
  }
  boost::shared_ptr<CompleteAbinitioMover> abinitio;
  if (stage_name == "complete") {
    score_fxn_name = "score3";
    core::scoring::ScoreFunctionOP score_4 = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score3").c_str());
    abinitio = boost::shared_ptr<CompleteAbinitioMover>(new CompleteAbinitioMover(score_4, de->frag_opt.frag_set_, de->frag_opt.frag_set_large));
  }

  //core::scoring::ScoreFunctionOP score = OptionsMapInitializer::CreateScoreFunctionStrategy::get()->build(std::string(score_fxn_name));
  //de->frag_opt.scorefxn = score;

  core::scoring::ScoreFunctionOP score = core::scoring::ScoreFunctionFactory::create_score_function(std::string(score_fxn_name).c_str());
  de->frag_opt.scorefxn = score;
  de->frag_opt.stage_name = stage_name;
  //  de->frag_mover = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::hybrid_mover, de->frag_opt);

  de->frag_mover = de->initialize_fragment_insertion_strategy(frag_type);

  if (stage_name == "complete") {
#if(USE_CRYO_EM)
    scfxn_ind_with_frags = boost::shared_ptr<FitFunction>( new PoseDensityFragmentFunction( de->pose_, score, de->ss, de->frag_mover));
    scfxn_ind_simple = boost::shared_ptr<FitFunction>( new PoseDensityFunction(de->pose_, score, de->ss, abinitio));
#else
    scfxn_ind_with_frags = boost::shared_ptr<FitFunction>( new PoseFragmentFunction(de->pose_, score, de->ss, abinitio));
    scfxn_ind_simple = boost::shared_ptr<FitFunction>( new PoseScoreFunction(de->pose_, score, de->ss, abinitio));
#endif
   } else {
#if(USE_CRYO_EM)
    scfxn_ind_with_frags = OptionsMapInitializer::FitnessFunctionInitializer(de->frag_opt, de->frag_mover, true).set_up();
     scfxn_ind_simple = OptionsMapInitializer::FitnessFunctionInitializer(de->frag_opt, de->frag_mover, true).set_up();
     //   scfxn_ind_with_frags = boost::shared_ptr<FitFunction>( new PoseDensityFragmentFunction(de->pose_, score, de->ss, de->frag_mover));
     //scfxn_ind_simple = boost::shared_ptr<FitFunction>( new PoseDensityFunction(de->pose_, score, de->ss, de->frag_mover));
#else
    scfxn_ind_with_frags = boost::shared_ptr<FitFunction>( new PoseFragmentFunction(de->pose_, score, de->ss, de->frag_mover));
    scfxn_ind_simple = boost::shared_ptr<FitFunction>( new PoseScoreFunction(de->pose_, score, de->ss, de->frag_mover));
#endif


  }


}



#endif /* UTILSMPI_H */
