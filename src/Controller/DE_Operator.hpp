#ifndef DEOPERATOR_HPP
#define DEOPERATOR_HPP
#include <core/pose/Pose.hh>
#include <boost/tokenizer.hpp>
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

#include "../Extra/Utils.hpp"
#include "../Algorithm/DifferentialEvolutionMover.hpp"
#include "../Movers/PoseFunction.hpp"
#include "../Movers/InitPopulation.hpp"
#include "../Movers/SecondStage.hpp"
#include "../Movers/FragInsertionMover.hpp"
#include "../Movers/CalculateRmsdDistancePopulation.hpp"
#include "../Movers/LocalSearchIndividualMover.hpp"
#include "FragmentInserterBuilder.hpp"
#include "OptionsMapInitializer.hpp"

class DE_Operator
{
public:
  std::map<std::string, std::string> score_per_stage;
  std::map<std::string, int> gmax_per_stage;
  std::map<std::string, Protinfo> prot_selection;
  std::map<std::string, OptionsMapInitializer::distances_enum> distances_map;
  std::map<std::string, OptionsMapInitializer::protocol_name_enum> protocol_name_map;
  std::map<std::string, OptionsMapInitializer::init_popul_strategy_enum> init_popul_strategy_map;

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

  CalculateDistancePopulationPtr use_distances_strategy(std::string option);

  void init_files(Protinfo prot_info);

  boost::shared_ptr<MoverDE> prepare_stage(std::string stage_name);

  void update_current_population_to_stage_score();

  void initialize_local_search_to_apply_at_population(std::string option);

  boost::shared_ptr<FragInsertionMover> initialize_fragment_insertion_strategy(std::string option );

  void initialize_population_stage1();

  boost::shared_ptr<MoverDE> init_differential_evolution_protocol(std::string option_protocol);

  void print_best_pose();

  void print_final_population(std::vector<Individual> popul );

  void apply_differential_evolution_with_stage(std::string stage_name);

  void init_two_stages_mover();

  void init_defaults();

  void init_available_stages();

  core::pose::PoseOP get_native_pose();

  virtual void run();

  void run_complete_abinitio();

public:
  int pose_size;
  int de_len;
  int start_res;
  double best_score;

  FragmentInserterBuilder frag_inserter_builder;

  std::string ss;
  core::pose::PoseOP pose_, best_pose_, native_pose_;
  core::scoring::ScoreFunctionOP scorefxn;
  core::pose::Pose ipose;
  protocols::moves::MonteCarloOP mc_;
  boost::shared_ptr<SecondStage > sstage;
  FragInsertionMoverPtr frag_mover;
  boost::shared_ptr<InitStagesMover> two_stages_mover;

  std::vector<Individual> current_population;
  FitFunctionPtr ffxn;
  PrintBestIndividualPtr print_best;
  CalculateDistancePopulationPtr calculate_distances_popul;
  InitPopulationPtr init_popul;
  boost::shared_ptr<LocalSearchIndividualMover> local_search;
  boost::property_tree::ptree app_options;
  std::vector<std::string> vec_stages;
};


#endif // DEoperator_hpp

