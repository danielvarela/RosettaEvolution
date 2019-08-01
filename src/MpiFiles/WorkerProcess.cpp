#include "WorkerProcess.hpp"

#include <boost/tokenizer.hpp>
#include "../Controller/DE_Operator.hpp"
#include "../Controller/StageBuilder.hpp"
#include <core/pose/copydofs/util.hh>
#include <core/pose/copydofs/CopyDofs.hh>
#include <core/pose/copydofs/CopyDofsInfo.hh>
#include <protocols/simple_moves/CopyDofMover.hh>
#include <core/pose/copydofs/util.hh>
#include "../Movers/PoseFunction.hpp"


WorkerProcess::ScoreStrategy::ScoreStrategy(   boost::shared_ptr<FitFunction> scfxn_ind_with_frags_1_in,  boost::shared_ptr<FitFunction> scfxn_ind_simple_1_in,
					       boost::shared_ptr<FitFunction> scfxn_ind_with_frags_2_in, boost::shared_ptr<FitFunction>  scfxn_ind_simple_2_in,
					       boost::shared_ptr<FitFunction> scfxn_ind_with_frags_3_in, boost::shared_ptr<FitFunction>  scfxn_ind_simple_3_in,
					       boost::shared_ptr<FitFunction> scfxn_ind_with_frags_complete_in, boost::shared_ptr<FitFunction>  scfxn_ind_simple_complete_in,
					       int mode_in, std::string stage_name_in
					       ) {
  scfxn_ind_with_frags_1 = scfxn_ind_with_frags_1_in;
  scfxn_ind_simple_1 = scfxn_ind_simple_1_in;
  scfxn_ind_with_frags_2 = scfxn_ind_with_frags_2_in;
  scfxn_ind_simple_2 = scfxn_ind_simple_2_in;
  scfxn_ind_with_frags_3 = scfxn_ind_with_frags_3_in;
  scfxn_ind_simple_3 = scfxn_ind_simple_3_in;
  scfxn_ind_with_frags_complete = scfxn_ind_with_frags_complete_in;
  scfxn_ind_simple_complete = scfxn_ind_simple_complete_in;
 #if(USE_CRYO_EM)
  
  scfxn_density= boost::dynamic_pointer_cast<PoseDensityFragmentFunction>(scfxn_ind_with_frags_3);

  boost::shared_ptr<PoseFragmentFunction> aux_ffxn = boost::dynamic_pointer_cast<PoseFragmentFunction>(scfxn_ind_with_frags_3);
 
  scfxn_density_simple = boost::shared_ptr<PoseDensityFunction>( new PoseDensityFunction(aux_ffxn->pose_, aux_ffxn->scorefxn, aux_ffxn->ss, aux_ffxn->frag_mover));
#endif
  mode = mode_in;
  stage_name = stage_name_in;
}

void
WorkerProcess::ScoreStrategy::apply(Individual& ind) {
  if (mode == 1) {
    double score = 0;
#if(USE_CRYO_EM)
    if (stage_name.compare("stage4")) score = scfxn_density_simple->score(ind);
#else
        if (stage_name.compare("stage2") == 0 ) score = scfxn_ind_simple_1->score(ind);
        if (stage_name.compare("stage3") == 0 ) {

          score = scfxn_ind_simple_2->score(ind);
            }
        if (stage_name.compare("stage4") == 0 ) {
          score = scfxn_ind_simple_3->score(ind);
        }
        //if (stage_name.compare("complete") == 0 ) score = scfxn_ind_simple_complete->score(ind);
        //std::cout << stage_name << " " << score << std::endl;
#endif
    
    ind.score = score;
  } else {
    if (mode == 2) {
      double score = 0;
#if(USE_CRYO_EM)
      score = scfxn_density->score(ind);
      #else
      if (stage_name == "stage2") score = scfxn_ind_with_frags_1->score(ind);
      if (stage_name == "stage3") {
        score = scfxn_ind_with_frags_2->score(ind);
      }
      if (stage_name == "stage4") score = scfxn_ind_with_frags_3->score(ind);
      //if (stage_name == "complete") score = scfxn_ind_with_frags_complete->score(ind);
#endif
      
      ind.score = score;
    }
  }
}



WorkerProcess::WorkerProcess(boost::shared_ptr<DE_Operator> app_operator) {
#if(USE_CRYO_EM)
  std::string frag_type = "hybrid_mover";
#else
  std::string frag_type = "stage_rosetta_mover";
#endif
  //std::cout << frag_type << std::endl;
  init_functions("stage2" ,frag_type, app_operator, scfxn_ind_with_frags_1, scfxn_ind_simple_1 );
  init_functions("stage3" ,frag_type, app_operator, scfxn_ind_with_frags_2, scfxn_ind_simple_2 );
  init_functions("stage4" ,frag_type, app_operator, scfxn_ind_with_frags_3, scfxn_ind_simple_3 );
  init_functions("complete" ,frag_type, app_operator, scfxn_ind_with_frags_complete, scfxn_ind_simple_complete );
#if(USE_CRYO_EM)
  fa_score = core::scoring::ScoreFunctionFactory::create_score_function("beta_cart");
  fa_score->set_weight( core::scoring::elec_dens_fast, 30.0 );
  to_cen = boost::shared_ptr<protocols::simple_moves::SwitchResidueTypeSetMover>( new protocols::simple_moves::SwitchResidueTypeSetMover("centroid"));

  to_all = boost::shared_ptr<protocols::simple_moves::SwitchResidueTypeSetMover>( new protocols::simple_moves::SwitchResidueTypeSetMover("fa_standard"));
  dockindens = protocols::electron_density::SetupForDensityScoringMoverOP( new protocols::electron_density::SetupForDensityScoringMover );
  my_relax = boost::shared_ptr<protocols::relax::FastRelax>( new protocols::relax::FastRelax(fa_score, 2 ));
  core::kinematics::MoveMapOP mm;
  mm = core::kinematics::MoveMapOP( new core::kinematics::MoveMap);
  mm->set_bb ( false );
  mm->set_chi ( true );
  mm->set_jump ( false );
  my_relax->set_movemap(mm);
#endif
 
}

void
WorkerProcess::run() {
  bool stop = true;
  while (stop) {
    int size;
    world.recv(0, 0, size);
    if (size == STOP_TAG) {
      break;
    }
    int mode;
    world.recv(0, 0, mode);
    std::string score_func;
    world.recv(0, 0, score_func);
    std::vector<IndMPI> popul(size);
    world.recv(0, 0, &popul.front(), popul.size() );
    if (mode == 66) {
      boost::shared_ptr<PoseFragmentFunction> pfunc = boost::dynamic_pointer_cast<PoseFragmentFunction>(scfxn_ind_with_frags_1);
      boost::shared_ptr<PoseDensityFunction> dens_func = boost::shared_ptr<PoseDensityFunction>( new PoseDensityFunction(pfunc->native, pfunc->scorefxn, pfunc->ss, NULL));

      for (int i = 0; i < popul.size(); i++) {
	core::pose::PoseOP my_pdb = pfunc->native->clone();
	//pfunc->fill_pose(my_pdb, popul[i].ind, popul[i].ind.ss);
	dockindens->apply(*my_pdb);
	dens_func->fill_pose(my_pdb, popul[i].ind, popul[i].ind.ss);

	to_all->apply(*my_pdb);
	my_relax->apply(*my_pdb);

	double result = (*fa_score)(*my_pdb);
	to_cen->apply(*my_pdb);
        dockindens->apply(*my_pdb);
 
	std::map<core::Size, core::Size > res_map;
	for (int j = 1; j <= my_pdb->total_residue() - 1; j++) {
	  res_map[j] = j;
	}
	protocols::simple_moves::CopyDofMover copy_dofs_in( *my_pdb, res_map);
	copy_dofs_in.apply( *my_pdb);
	core::pose::copydofs::CopyDofsInfo dofs_info_container = copy_dofs_in.copy_dofs_info( *my_pdb );
	utility::vector1< std::pair< core::id::DOF_ID, core::Real > > vector_info =  dofs_info_container.dofs_info();

	std::vector<double> dofs_vars;
	for ( auto const & elem : vector_info ) {
	  core::Real dof_value = elem.second; // Index in the little "chunk" or "scratch" pose
	  dofs_vars.push_back(dof_value);
	}
	popul[i].ind.vars = dofs_vars;

	// old version
	//pfunc->pose_to_ind(my_pdb, popul[i].ind);
	popul[i].ind.score = result;


      }

    } else {
      ScoreStrategy score_ind( scfxn_ind_with_frags_1,  scfxn_ind_simple_1,
			     scfxn_ind_with_frags_2,  scfxn_ind_simple_2,
			     scfxn_ind_with_frags_3,  scfxn_ind_simple_3,
			     scfxn_ind_with_frags_complete,  scfxn_ind_simple_complete,
			     mode, score_func);
      score_ind.stage_name = score_func;
      for (int i = 0; i < popul.size(); i++) {
	score_ind.apply(popul[i].ind);
      }

    }

    world.send(0, 0, &popul.front(), popul.size() ) ;
  }
}



WorkerProcessEvaluateAndNearest::WorkerProcessEvaluateAndNearest(boost::shared_ptr<DE_Operator> app_operator) : WorkerProcess(app_operator) {
  StageBuilder stage_builder(app_operator->app_options, app_operator->scorefxn,  app_operator->frag_opt);
  calculate_distances_popul = CalculateRmsdDistancePopulation(app_operator->native_pose_, app_operator->ffxn, app_operator->ss, app_operator->scorefxn, app_operator->fit_radius ).use_distances_strategy( app_operator->app_options.get<std::string>("Protocol.distance_strategy") );

  // stage_builder.initializer_frag_mover();
  // stage_builder.initializer_frag_popul();
  // stage_builder.initializer_print_best();
  // stage_builder.initializer_distances_calculator();
  // calculate_distances_popul = stage_builder.calculate_distances_popul;

}


void
WorkerProcessEvaluateAndNearest::run() {
  bool stop = true;
  int TAG_ERROR =  150;
  while (stop) {
    int popul_size;
    world.recv(0, 0, popul_size);
    if (popul_size == STOP_TAG) {
      break;
    }
    // std::cout << "popul size " << popul_size << std::endl;
    std::vector<IndMPI> previous_popul(popul_size);
    world.recv(0, 0, &previous_popul.front(), previous_popul.size() );
    std::vector<Individual> right_popul;
    for (int i = 0; i < previous_popul.size(); i++) {
      right_popul.push_back(previous_popul[i].ind);
    }

    int size;
    world.recv(0, 0, size);
    int mode;
    world.recv(0, 0, mode);
    std::string score_func;
    world.recv(0, 0, score_func);

    // acaba de recibir los datos en comun
    IndMPI ind_received;
    world.recv(0, 0, ind_received );

    ScoreStrategy score_ind( scfxn_ind_with_frags_1,  scfxn_ind_simple_1,
			     scfxn_ind_with_frags_2,  scfxn_ind_simple_2,
			     scfxn_ind_with_frags_3,  scfxn_ind_simple_3,
			     scfxn_ind_with_frags_complete,  scfxn_ind_simple_complete,
			     mode, score_func);


    while (ind_received.ind.vars[0] != TAG_ERROR) {
      // score
      score_ind.apply(ind_received.ind);
      // nearest
      int nearest_ind = calculate_distances_popul->find_nearest(ind_received.ind, right_popul );
      ind_received.nearest = nearest_ind;
      // send back the information as a IndMPI (including nearest index for each ind)
      world.send(0, 0, ind_received );
      world.recv(0, 0, ind_received );
    }
  }
}


WorkerProcessScatterGather::WorkerProcessScatterGather(boost::shared_ptr<DE_Operator> app_operator) : WorkerProcess(app_operator) {
  StageBuilder stage_builder(app_operator->app_options, app_operator->scorefxn,  app_operator->frag_opt);
  calculate_distances_popul = CalculateRmsdDistancePopulation(app_operator->native_pose_, app_operator->ffxn, app_operator->ss, app_operator->scorefxn, app_operator->fit_radius ).use_distances_strategy( app_operator->app_options.get<std::string>("Protocol.distance_strategy") );
}

void
WorkerProcessScatterGather::run() {
  bool stop = true;
  int TAG_ERROR =  150;
  while (stop) {
    int popul_size;
    world.recv(0, 0, popul_size);
    if (popul_size == STOP_TAG) {
      exit(0);
      break;
    }
    // std::cout << "popul size " << popul_size << std::endl;
    std::vector<IndMPI> previous_popul(popul_size);
    world.recv(0, 0, &previous_popul.front(), previous_popul.size() );
    std::vector<Individual> right_popul;
    for (int i = 0; i < previous_popul.size(); i++) {
      right_popul.push_back(previous_popul[i].ind);
    }

    int size;
    world.recv(0, 0, size);
    int mode;
    world.recv(0, 0, mode);
    std::string score_func;
    world.recv(0, 0, score_func);
    // acaba de recibir los datos comunes

    ScoreStrategy score_ind( scfxn_ind_with_frags_1,  scfxn_ind_simple_1,
			     scfxn_ind_with_frags_2,  scfxn_ind_simple_2,
			     scfxn_ind_with_frags_3,  scfxn_ind_simple_3,
			     scfxn_ind_with_frags_complete,  scfxn_ind_simple_complete,
			     mode, score_func);

    std::vector<std::vector<IndMPI> > vector_parts;
    std::vector<IndMPI> sub_popul;
    scatter(world, vector_parts, sub_popul, 0);
    for (int i = 0; i < sub_popul.size(); i++) {
      score_ind.apply(sub_popul[i].ind);
      int nearest_ind = calculate_distances_popul->find_nearest(sub_popul[i].ind, right_popul );
      sub_popul[i].nearest = nearest_ind;
    }
    gather(world, sub_popul, vector_parts, 0);

  }
}



/*


   timestamp_t init_ScoreCalculation = get_timestamp();

    for (int i = 0; i < trial_popul.size(); i++) {
      score_ind.apply(trial_popul[i].ind);
    }
    timestamp_t end_ScoreCalculation = get_timestamp();
    double secs_ScoreCalculation = (end_ScoreCalculation - init_ScoreCalculation) / 1000000.0L;
    std::cout << "NearCalcualtion time " << secs_ScoreCalculation << std::endl;


    // find nearest of each trial_popul i
    timestamp_t init_NearCalculation = get_timestamp();

    int nearest_ind;
    for (int i = 0; i < trial_popul.size(); i++) {
      nearest_ind = calculate_distances_popul->find_nearest(trial_popul[i].ind, right_popul );
      // nearest_ind = 0;
      trial_popul[i].nearest = nearest_ind;
    }
    timestamp_t end_NearCalculation = get_timestamp();
    double secs_NearCalculation = (end_NearCalculation - init_NearCalculation) / 1000000.0L;
    std::cout << "NearCalcualtion time " << secs_NearCalculation << std::endl;



 */
