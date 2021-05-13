
#include "MasterRosettaCalculator.hpp"
#include <vector>
#include "../Controller/DE_Operator.hpp"
#include <core/pose/copydofs/util.hh>
#include <core/pose/copydofs/CopyDofs.hh>
#include <core/pose/copydofs/CopyDofsInfo.hh>
#include <protocols/simple_moves/CopyDofMover.hh>
#include <core/pose/copydofs/util.hh>




MasterRosettaCalculator::MasterRosettaCalculator(boost::shared_ptr<FitFunction> scfxn_ind_with_frags_input)  {
  scfxn_ind_with_frags = scfxn_ind_with_frags_input;
  boost::shared_ptr<PoseFragmentFunction> aux_ffxn = boost::dynamic_pointer_cast<PoseFragmentFunction>(scfxn_ind_with_frags);
  

  simple_ffxn = boost::shared_ptr<PoseScoreFunction>( new PoseScoreFunction(aux_ffxn->pose_, aux_ffxn->scorefxn, aux_ffxn->ss, aux_ffxn->frag_mover));
#if(USE_CRYO_EM)
  simple_ffxn = boost::shared_ptr<PoseScoreFunction>( new PoseDensityFunction(aux_ffxn->pose_, aux_ffxn->scorefxn, aux_ffxn->ss, aux_ffxn->frag_mover));
  density_func = boost::shared_ptr<PoseDensityFunction>( new PoseDensityFunction(aux_ffxn->pose_, aux_ffxn->scorefxn, aux_ffxn->ss, aux_ffxn->frag_mover));
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
#else
  simple_ffxn = boost::shared_ptr<PoseScoreFunction>( new PoseScoreFunction(aux_ffxn->pose_, aux_ffxn->scorefxn, aux_ffxn->ss, aux_ffxn->frag_mover));
#endif
}

std::vector<std::vector<IndMPI> >
MasterRosettaCalculator::SplitVector(const std::vector<Individual>& vec, size_t n)
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
MasterRosettaCalculator::run(std::vector<Individual> num_popul, int mode) {
  std::vector<boost::mpi::request> reqs;
  // desde aqui es la funcion calcular popul
  //std::string score_name = simple_ffxn->name();
  std::string score_name = stage_name;
  std::string score_func = stage_name;
  // if (score_name == "complete") {
  //   score_func = "complete";
  // } else {
  //   std::cout << "SCORE NAME!! " << score_name << std::endl;
  //   //score_func = std::string(score_name.substr(score_name.rfind(" ") + 1));
  //   std::cout << "Score func " << score_func << std::endl;
  // }
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
#if(USE_CRYO_EM)
      density_func->score(master_vector_for_work[i].ind);
#else
      simple_ffxn->score(master_vector_for_work[i].ind);
#endif
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


std::vector<Individual>
MasterRosettaCalculator::dump(std::vector<Individual> num_popul) {
  std::vector<boost::mpi::request> reqs;
  int mode = 66;
  // desde aqui es la funcion calcular popul
  std::string score_name = scfxn_ind_with_frags->name();
  std::string score_func;
  score_func = "complete";
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
  boost::shared_ptr<PoseFragmentFunction> pfunc = boost::dynamic_pointer_cast<PoseFragmentFunction>(scfxn_ind_with_frags);
  core::pose::PoseOP native = boost::dynamic_pointer_cast<PoseFragmentFunction>(scfxn_ind_with_frags)->native->clone();
  boost::shared_ptr<PoseDensityFunction> dens_func = boost::shared_ptr<PoseDensityFunction>( new PoseDensityFunction(pfunc->native, pfunc->scorefxn, pfunc->ss, NULL));


  std::cout << "llega a dump popul " << std::endl;
  for (int i = 0; i < master_vector_for_work.size(); i++) {
    core::pose::PoseOP my_pdb = native->clone();
    //pfunc->fill_pose(my_pdb, master_vector_for_work[i].ind, master_vector_for_work[i].ind.ss);
    dockindens->apply(*my_pdb);
    dens_func->fill_pose(my_pdb, master_vector_for_work[i].ind, master_vector_for_work[i].ind.ss);
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
    master_vector_for_work[i].ind.vars = dofs_vars;

    //pfunc->pose_to_ind(my_pdb, master_vector_for_work[i].ind);
    master_vector_for_work[i].ind.score = result;
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




std::vector<IndMPI>
MasterEvaluateAndNearestCalculator::run(  std::vector<Individual> popul, std::vector<Individual> trial_popul, int mode) {
  std::vector<boost::mpi::request> reqs;
  // desde aqui es la funcion calcular popul
  std::string score_name = scfxn_ind_with_frags->name();
  std::string score_func;
  if (score_name == "complete") {
    score_func = "complete";
  } else {
    score_func = std::string(score_name.substr(score_name.rfind(" ") + 1));
  }

  std::vector<IndMPI> mpi_popul;
  for (int i = 0; i < popul.size(); i++) {
    mpi_popul.push_back(IndMPI(i, popul[i]));
  }
  IndMPI dummy_individual = mpi_popul[0];
  int TAG_ERROR =  150;
  dummy_individual.ind.vars[0] = TAG_ERROR;

  int array_size = trial_popul.size() / (world.size());
  int mpi_popul_size =  mpi_popul.size();
  //int mode = 2; // for frag insertion
  std::vector<IndMPI> master_vector_for_work;
  std::vector<std::vector<IndMPI> > vector_parts = SplitVector(trial_popul, world.size());
  //  master_vector_for_work = vector_parts[0];

  //  envio data comun
  for (int i = 1; i < world.size(); i++) {
    world.send(i, 0, mpi_popul_size);
    world.send(i, 0, &mpi_popul.front(), mpi_popul.size() );
    world.send(i, 0, array_size);
    world.send(i, 0, mode);
    world.send(i, 0, score_func);
  }

  std::vector<boost::mpi::request> reqs_receive(world.size() - 1);
  std::vector<IndMPI> mpi_result_popul(popul.size());
  int num_processed_inds = 0;
  int c_pending = 0;
  for (int i = 1; i < world.size(); i++) {
    world.send(i, 0, mpi_popul[num_processed_inds] );
    c_pending++;
    reqs_receive[i - 1] = world.irecv(i, 0, mpi_result_popul[num_processed_inds] );
    num_processed_inds++;
  }

  while (c_pending != 0) {
    std::pair<boost::mpi::status, std::vector<boost::mpi::request>::iterator > return_value = mpi::wait_any( std::begin(reqs_receive) , std::end(reqs_receive) );
    int process_that_sended_data = return_value.first.source();
    c_pending--;
    if (num_processed_inds < popul.size()) {
      world.send(process_that_sended_data, 0, mpi_popul[num_processed_inds] );
      c_pending++;
      *return_value.second = world.irecv( process_that_sended_data, 0, mpi_result_popul[num_processed_inds] );
      num_processed_inds++;
    } else {
      world.send(process_that_sended_data, 0, dummy_individual ) ;
    }
  }

  std::vector<IndMPI> my_new_result(mpi_result_popul.size());
  for (int i = 0; i < mpi_result_popul.size(); i++) {
    my_new_result[mpi_result_popul[i].index] = mpi_result_popul[i];
  }

  return my_new_result;
}




std::vector<IndMPI>
MasterScatterGather::run(  std::vector<Individual> popul, std::vector<Individual> trial_popul, int mode) {
  std::vector<boost::mpi::request> reqs;
  // desde aqui es la funcion calcular popul
  std::string score_name = scfxn_ind_with_frags->name();
  std::string score_func;
  if (score_name == "complete") {
    score_func = "complete";
  } else {
    score_func = std::string(score_name.substr(score_name.rfind(" ") + 1));
  }

  std::vector<IndMPI> mpi_popul;
  for (int i = 0; i < popul.size(); i++) {
    mpi_popul.push_back(IndMPI(i, popul[i]));
  }
  IndMPI dummy_individual = mpi_popul[0];
  int TAG_ERROR =  150;
  dummy_individual.ind.vars[0] = TAG_ERROR;

  int array_size = trial_popul.size() / (world.size());
  int mpi_popul_size =  mpi_popul.size();
  //int mode = 2; // for frag insertion
  std::vector<IndMPI> master_vector_for_work;
  std::vector<std::vector<IndMPI> > vector_parts = SplitVector(trial_popul, world.size());
  //  master_vector_for_work = vector_parts[0];

  //  envio data comun
  for (int i = 1; i < world.size(); i++) {
    world.send(i, 0, mpi_popul_size);
    world.send(i, 0, &mpi_popul.front(), mpi_popul.size() );
    world.send(i, 0, array_size);
    world.send(i, 0, mode);
    world.send(i, 0, score_func);
  }

  std::vector<IndMPI> sub_popul;
  scatter(world, vector_parts, sub_popul, 0);
  for (int i = 0; i < sub_popul.size(); i++) {
    if (mode == 1) {
      simple_ffxn->score(sub_popul[i].ind);
    } else {
      scfxn_ind_with_frags->score(sub_popul[i].ind);
    }
    int nearest_ind = calculate_distances_popul->find_nearest(sub_popul[i].ind, popul );
    sub_popul[i].nearest = nearest_ind;
  }
  gather(world, sub_popul, vector_parts, 0);

  std::vector<IndMPI> mpi_result_popul(popul.size());
  for (int i = 0; i < vector_parts.size(); i++) {
    for (std::vector<IndMPI>::iterator it = vector_parts[i].begin(); it != vector_parts[i].end(); ++it) {
      mpi_result_popul[it->index] = *it;
    }
  }

  return mpi_result_popul;
}
