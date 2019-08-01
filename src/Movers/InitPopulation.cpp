
#include "InitPopulation.hpp"

void
InitPopulation::apply(std::vector<Individual>& popul, int NP, int ind_size) {

  for (int i = 0; i < NP; ++i) {
    Individual ind(ind_size);

    for (int j = 0; j < ind_size; ++j) {
      if (URAND() < 0.5) {
	ind.vars[j] = URAND();
      } else {
	ind.vars[j] = URAND() * -1;
      }
    }

    ind.vars[0] = 0;
    ind.vars[ind.vars.size() - 1] = 0;

    score(ind);

    popul.push_back(ind);
  }
}



void
PosePopulation::apply(std::vector<Individual>& popul, int NP, int ind_size) {
  //    std::cout << "ind_size " <<  ind_size << std::endl;
  popul.resize(0);
  for (int i = 0; i < NP; ++i) {
    Individual ind(ind_size);
    disturb_individual(ind, ind_size);
    // ind.vars[0] = 0;
    // ind.vars[ind.vars.size() - 1] = 0;
    score(ind);
    popul.push_back(ind);
  }
}



void
PosePopulation::disturb_individual(Individual &ind, int ind_size) {
  int nres = scfxn->start_res();
  //dist_degrees = 5;
  // for (int j = 0; j < ind_size - 2; j = j+2) {
  //   //taking into account the previous obtained angles
  //   if (URAND() < 0.5) {
  //     if (URAND() > 0.5) {
  // 	ind.vars[j] =  pose_->phi(nres) + (URAND() * dist_degrees) ;
  //     } else {
  // 	ind.vars[j] = pose_->phi(nres) - (URAND() * dist_degrees);
  //     }
  //     if (ind.vars[j] > 180.0) {
  // 	ind.vars[j] -= std::abs(ind.vars[j] - 180.0);
  //     }
  //     if (ind.vars[j] < -180.0) {
  // 	ind.vars[j] += std::abs(ind.vars[j] - 180.0);
  //     }
  //   } else {
  //     ind.vars[j] = pose_->phi(nres);
  //   }
  //   ind.vars[j] = scale(ind.vars[j]);
  //   if (URAND() < 0.5) {
  //     if (URAND() > 0.5) {
  // 	ind.vars[j + 1] = pose_->psi(nres) + (URAND() * dist_degrees) ;
  //     } else {
  // 	ind.vars[j + 1] = pose_->psi(nres) - (URAND() * dist_degrees);
  //     }

  //     if (ind.vars[j + 1] > 180.0) {
  // 	ind.vars[j + 1] -= std::abs(ind.vars[j + 1] - 180.0);
  //     }
  //     if (ind.vars[j + 1] < -180.0) {
  // 	ind.vars[j + 1] += std::abs(ind.vars[j + 1] - 180.0);
  //     }
  //   }else {
  //     ind.vars[j + 1] = pose_->psi(nres);
  //   }
  //   ind.vars[j+1] = scale(ind.vars[j+1]);
  //   nres++;
  // }


  core::pose::PoseOP inner_pose = pose_->clone();
  for (int i = 1; i <= pose_->total_residue() - 1; i++) {
     if (URAND() > 0.5) {
       inner_pose->set_phi(i , pose_->phi(i) + (URAND() * dist_degrees) ) ;
      } else {
       inner_pose->set_phi(i , pose_->phi(i) - (URAND() * dist_degrees) ) ;
      }
      if (URAND() > 0.5) {
       inner_pose->set_psi(i , pose_->psi(i) + (URAND() * dist_degrees) ) ;
      } else {
       inner_pose->set_psi(i , pose_->psi(i) - (URAND() * dist_degrees) ) ;
      }
  }

  //inner_pose->dump_pdb("example_popul_pose.pdb");
  
    std::map<core::Size, core::Size > res_map;
    for (int i = 1; i <= inner_pose->total_residue() - 1; i++) {
      res_map[i] = i;
    }
    protocols::simple_moves::CopyDofMover copy_dofs_in( *inner_pose, res_map);
    copy_dofs_in.apply( *inner_pose );
    core::pose::copydofs::CopyDofsInfo dofs_info_container = copy_dofs_in.copy_dofs_info( *inner_pose );
    utility::vector1< std::pair< core::id::DOF_ID, core::Real > > vector_info =  dofs_info_container.dofs_info();
    std::vector<double> dofs_vars;
    for ( auto const & elem : vector_info ) {
      core::Real dof_value = elem.second; // Index in the little "chunk" or "scratch" pose
      dofs_vars.push_back(dof_value);
      //if (dofs_vars.size() < 10) {
      //	std::cout << dof_value << " , ";
      //}
    }
    // std::cout << std::endl;
    ind.vars = dofs_vars;

  ind.omega.resize(0);
  for (int i = 1; i < pose_->total_residue(); i++) {
    ind.omega.push_back(pose_->omega(i));
  }
  ind.ss = pose_->secstruct();
  ind.score = 1000;
};

void
PosePopulation::secstruct_individual(Individual &ind, int ind_size) {

  std::string ss_struct = pose_->secstruct();
  int nres = scfxn->start_res();
  for (int j = 0; j < ind_size - 2; j = j+2) {
    // secstruct phi and psi
    if (ss_struct[nres - 1] == 'H') {
      if (URAND() > 0.5) {
	ind.vars[j] = scale(-64 + (URAND() * 7));
      } else {
	ind.vars[j] = scale(-64 - (URAND() * 7));
      }

      if (URAND() > 0.5) {
	ind.vars[j + 1] = scale(-41 + (URAND() * 7)  );
      } else {
	ind.vars[j + 1] = scale(-41 - (URAND() * 7)  );
      }

    } else {

      if (URAND() < 0.5) {
	ind.vars[j] = URAND();
      } else {
	ind.vars[j] = URAND() * -1;
      }

    }


    nres++;
  }

}

void
PosePopulation::rama_individual(Individual &ind, int ind_size) {
  core::scoring::Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();

  std::string ss_struct = pose_->secstruct();

  double rama_phi, rama_psi;

  int nres = scfxn->start_res();
  for (int j = 0; j < ind_size - 2; j = j+2) {

    // rama phi and psi
    rama.random_phipsi_from_rama(pose_->aa(nres), rama_phi, rama_psi);
    ind.vars[j] = scale_rama(rama_phi);
    ind.vars[j + 1] = scale_rama(rama_psi);
    //	std::coutmo << "nres " << nres << " j " << j << " | " << rama_phi << " " << rama_psi << std::endl;
    nres++;
  }
}

double
PosePopulation::scale_rama( double old_value) {
  double old_min = 0.0;
  double old_max = 360.0;
  double new_min = -1*DE_LIMIT;
  double new_max = DE_LIMIT;

  return ( (old_value - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min;
}

double
PosePopulation::scale( double old_value) {
  double old_min = -180.0;
  double old_max = 180.0;
  double new_min = -1.0*DE_LIMIT;
  double new_max = DE_LIMIT;

  return ( (old_value - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min;
}

void
PoseTotalRandomInitPopulation::apply(std::vector<Individual>& popul, int NP, int ind_size) {

  boost::shared_ptr<PoseFragmentFunction> pfunc = boost::dynamic_pointer_cast<PoseFragmentFunction >(scfxn);
  for (int i = 0; i < NP; ++i) {
    Individual ind(ind_size);
    disturb_individual(ind, ind_size);
    pfunc->fill_pose(pose_, ind, ss);
    pfunc->pose_to_ind(pose_, ind );
    score(ind);
    popul.push_back(ind);
  }

}

void
PoseTotalRandomInitPopulation::disturb_individual(Individual &ind, int ind_size) {
  ind.vars.resize(ind_size);
  for (int j = 0; j < ind_size; ++j) {
    if (URAND() < 0.5) {
      ind.vars[j] = URAND();
    } else {
      ind.vars[j] = URAND() * -1;
    }
  }

  ind.vars[0] = 0;
  ind.vars[ind.vars.size() - 1] = 0;
  ind.omega.resize(0);
  for (int i = 1; i < pose_->total_residue(); i++) {
    ind.omega.push_back(pose_->omega(i));
  }
  ind.ss = pose_->secstruct();
  ind.score = 1000;
}

void
TwoStagesInitPopulation::apply(std::vector<Individual>& popul, int NP, int ind_size) {
  boost::shared_ptr<PoseScoreFunction> pfunc = boost::dynamic_pointer_cast<PoseScoreFunction >(scorefxn);
  for (int i = 0; i < NP; ++i) {
    Individual ind(ind_size);
    disturb_individual(ind, ind_size);
    pfunc->fill_pose(my_pose_, ind, ss);
    two_stages_mover->apply(*my_pose_);
    pfunc->pose_to_ind(my_pose_, ind );
    ind.vars[0] = 0;
    ind.vars[ind.vars.size() - 1] = 0;
    scorefxn->score(ind);
    popul.push_back(ind);
  }

}

void
TwoStagesInitPopulation::disturb_individual(Individual &ind, int ind_size) {
  ind.vars.resize(ind_size);
  for (int j = 0; j < ind_size; ++j) {
    if (URAND() < 0.5) {
      ind.vars[j] = URAND();
    } else {
      ind.vars[j] = URAND() * -1;
    }
  }

  ind.vars[0] = 0;
  ind.vars[ind.vars.size() - 1] = 0;
  ind.omega.resize(0);
  for (int i = 1; i < my_pose_->total_residue(); i++) {
    ind.omega.push_back(my_pose_->omega(i));
  }
  ind.ss = my_pose_->secstruct();
  ind.score = 1000;
}
