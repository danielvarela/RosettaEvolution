#ifndef DEOPERATOR_HPP
#define DEOPERATOR_HPP

#include <core/pose/Pose.hh>

#include "Utils.hpp"
#include "moves/DifferentialEvolutionMover.hpp"
#include "moves/PoseFunction.hpp"
#include "moves/InitPopulation.hpp"
#include "moves/SecondStage.hpp"
#include "moves/FragInsertionMover.hpp"

#include "moves/CalculateRmsdDistancePopulation.hpp"
#include "moves/LocalSearchIndividualMover.hpp"
#include <boost/tokenizer.hpp>

//#include "ReaderFragmentBuilder.hpp"
//#include "relax.hpp"

#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/elec/FA_ElecEnergy.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/lkball/LK_BallEnergy.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/nonlocal/SingleFragmentMover.hh>


#include <boost/property_tree/ptree.hpp>


#include <algorithm>
#include <cctype>
#include <locale>

// trim from start (in place)
static inline void ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
	return !std::isspace(ch);
      }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
	return !std::isspace(ch);
      }).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
  ltrim(s);
  rtrim(s);
}


class DE_Operator
{
public:
  std::map<std::string, std::string> score_per_stage;
  std::map<std::string, int> gmax_per_stage;
  std::map<std::string, Protinfo> prot_selection;

  enum distances_enum {
    rmsd, rmsd_native_diff, euclidean, euclidean_loop, euclidean_diff_abs, euclidean_partial_mario, euclidean_mario, rmsd_without_superp
  };
  std::map<std::string, distances_enum> distances_map;

  enum protocol_name_enum{
    Shared, HybridShared
  };
  std::map<std::string, protocol_name_enum> protocol_name_map;

  enum init_popul_strategy_enum{
    total_random, total_random_pose_based, random_pose_based, init_popul_with_stage
  };
  std::map<std::string, init_popul_strategy_enum> init_popul_strategy_map;

  FragInsertionStrategy::FragOptions frag_opt;
  FragInsertionMoverPtr local_search_frag_mover;
  double fit_radius;

  boost::property_tree::ptree get_app_options() {
    return app_options;
  }

  DE_Operator();

  DE_Operator(std::string prot_selected, boost::property_tree::ptree opt);

  DE_Operator(std::string prot_selected);

  void init_protein_archive();

  void init_setup();

  boost::shared_ptr<InitPopulation> initialize_init_popul_strategy(std::string init_popul_option );

  CalculateRmsdDistancePopulationPtr use_distances_strategy(std::string option);

  void init_files(Protinfo prot_info);

  boost::shared_ptr<MoverDE> prepare_stage(std::string stage_name);

  void update_current_population_to_stage_score();

  void initialize_population_stage1();

  boost::shared_ptr<MoverDE> init_differential_evolution_protocol();

  void print_best_pose();

  void apply_differential_evolution_with_stage(std::string stage_name);

  void init_two_stages_mover();

  void init_default_score();
 
  void init_available_stages();

  core::pose::PoseOP get_native_pose();

  virtual void run();

public:
  int pose_size;
  int de_len;
  int start_res;
  double best_score;
  std::string ss;
  core::pose::PoseOP pose_, best_pose_, native_pose_;
  core::scoring::ScoreFunctionOP scorefxn;
  core::pose::Pose ipose;
  protocols::moves::MonteCarloOP mc_;
  boost::shared_ptr<SecondStage > sstage;
  FragInsertionMoverPtr frag_mover;
  boost::shared_ptr<InitStagesMover> two_stages_mover;
  core::fragment::FragSetOP frag_set_, frag_set_large;

  std::vector<Individual> current_population;
  FitFunctionPtr ffxn;
  PrintBestIndividualPtr print_best;
  CalculateRmsdDistancePopulationPtr calculate_distances_popul;
  InitPopulationPtr init_popul;
  boost::shared_ptr<LocalSearchIndividualMover> local_search;
  boost::property_tree::ptree app_options;
  std::vector<std::string> vec_stages;
};


#endif // DEoperator_hpp

/*

Bakcup of my init setup of DE Mover, but its almost the same code as the function
                boost::shared_ptr<MoverDE> prepare_stage(stage_name);

    // FitFunctionPtr ffxn = boost::shared_ptr<FitFunction>( new PoseFunction(pose_, scorefxn, ss));
    // FitFunctionPtr ffxn = boost::shared_ptr<FitFunction>( new PoseFragmentFunction(pose_, scorefxn, ss, frag_mover));
    // PrintBestIndividualPtr print_best = PrintBestIndividualPtr( new PrintBestIndividual(pose_, ffxn, ss, scorefxn));
    // CalculateRmsdDistancePopulationPtr calculate_distances_popul = CalculateRmsdDistancePopulationPtr( new CalculateRmsdDistancePopulation(pose_, ffxn, ss, scorefxn));
    // InitPopulationPtr init_popul = boost::shared_ptr<InitPopulation>(new PosePopulation(ffxn, pose_, ss));

    // //    InitPopulationPtr init_popul = boost::shared_ptr<InitPopulation>(new TwoStagesInitPopulation(ffxn, pose_, two_stages_mover, ss));


 */


/*COMPLETE ABINITIO MOVER
    // int cnt = 0;
    // do {
    //   pose_ = native_->clone();
    // frag_mover = FragInsertionMoverPtr( new CompleteAbinitioMover(scorefxn, frag_set_, frag_set_large));
    // frag_mover->apply(*pose_);
    // std::cout << "FINAL SCORE " << cnt << " : " << (*scorefxn)(*pose_) << std::endl;
    // std::cout << "rmsd : " << core::scoring::CA_rmsd(*pose_, *native_) << std::endl;
    // cnt++;
    // } while (cnt < 100);
    //    sstage = boost::shared_ptr<SecondStage >( new SecondStage() );


 */
