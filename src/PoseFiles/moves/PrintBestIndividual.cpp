
#include "PrintBestIndividual.hpp"
#include <string>
#include <vector>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>


PrintBestIndividual::PrintBestIndividual(const core::pose::PoseOP& p , FitFunctionPtr sfxn,  std::string ss_in, core::scoring::ScoreFunctionOP pscore)  {
  pose_ = p->clone();
  uniq_id = static_cast<int>(std::abs(rand()));
  cnt = 0;
  scorefxn = sfxn;
  ss = ss_in;
  native_ = p->clone();
  pose_score = pscore;
}

std::string
PrintBestIndividual::print(Individual ind) {
 return print(ind, std::to_string(uniq_id) + "_" + std::to_string(cnt) );
}

std::string
PrintBestIndividual::print(Individual ind, std::string id) {
  boost::shared_ptr<PoseFunction> pfunc = boost::dynamic_pointer_cast<PoseFunction >(scorefxn);
  pfunc->fill_pose(pose_, ind, ss);
  std::string path = "/home/dvarela/Code/RosettaEvolution/output_pdbs/best_ind_" + id +"_.pdb";

#if(USE_CRYO_EM)
  protocols::electron_density::SetupForDensityScoringMoverOP dockindens;
  dockindens = protocols::electron_density::SetupForDensityScoringMoverOP( new protocols::electron_density::SetupForDensityScoringMover );
     dockindens->apply(*pose_);
   #endif

     pose_->dump_pdb(path.c_str());

  double rmsd_vs_previous_best = 0;
  if (cnt == 0) {
    rmsd_vs_previous_best = 0;
    previous_best = pose_->clone();
  } else {
    if (previous_best) {
      rmsd_vs_previous_best = core::scoring::CA_rmsd(*previous_best, *pose_);
    } else {
      rmsd_vs_previous_best = 0;
    }
    previous_best = pose_->clone();

  }
  path = path + " " + std::to_string(core::scoring::CA_rmsd(*pose_, *native_));
  path = path + " " + std::to_string(rmsd_vs_previous_best);
  cnt++;
  return path;
}


std::string
PrintBestIndividual::print_energies(Individual ind) {
  boost::shared_ptr<PoseFunction> pfunc = boost::dynamic_pointer_cast<PoseFunction >(scorefxn);
  pfunc->fill_pose(pose_, ind, ss);
  (*pose_score)(*pose_);
  //pose_score->show(std::cout, *pose_);
  core::scoring::EnergyMap energy_map = pose_->energies().total_energies();
  std::string energies_string;
  energies_string += "env " + std::to_string(energy_map[core::scoring::ScoreType::env]) + " ";
  energies_string += "pair " + std::to_string(energy_map[core::scoring::ScoreType::pair]) + " ";
  if (pose_score->get_weight(core::scoring::ScoreType::cbeta) != 0 ) {
    energies_string += "cbeta " + std::to_string(energy_map[core::scoring::ScoreType::cbeta]) + " ";
  }

  energies_string += "vdw " + std::to_string(energy_map[core::scoring::ScoreType::vdw]) + " ";
  if (pose_score->get_weight(core::scoring::ScoreType::rg) != 0 ) {
    energies_string += "rg " + std::to_string(energy_map[core::scoring::ScoreType::rg]) + " ";
  }

  if (pose_score->get_weight(core::scoring::ScoreType::cenpack) != 0 ) {
    energies_string += "cenpack " + std::to_string(energy_map[core::scoring::ScoreType::cenpack]) + " ";
  }
  energies_string += "hs_pair " + std::to_string(energy_map[core::scoring::ScoreType::hs_pair]) + " ";
  energies_string += "ss_pair " + std::to_string(energy_map[core::scoring::ScoreType::ss_pair]) + " ";
  if (pose_score->get_weight(core::scoring::ScoreType::rsigma) != 0 ) {
    energies_string += "rsigma " + std::to_string(energy_map[core::scoring::ScoreType::rsigma]) + " ";
  }

  energies_string += "sheet " + std::to_string(energy_map[core::scoring::ScoreType::sheet]) + " ";
  return energies_string;
}
