

#ifndef UTILS_HPP
#define UTILS_HPP


#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/pose/util.hh>

//Density
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/elec/FA_ElecEnergy.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/lkball/LK_BallEnergy.hh>
#include <core/fragment/FragmentIO.hh>




struct Protinfo
{
  std::string name;
  std::string pdb_file;
  std::string map_file;
  std::string ss_file;
  std::string frag_3;
  std::string frag_9;
};
#include <sys/types.h>
#include <sys/stat.h>


inline int dirExists(const char *path)
{
  struct stat info;

  if(stat( path, &info ) != 0)
    return 0;
  else if(info.st_mode & S_IFDIR)
    return 1;
  else
    return 0;
}
inline std::string read_ss2(std::string ss_file, core::pose::Pose& ipose) {

  core::fragment::SecondaryStructureOP ss_profile( new core::fragment::SecondaryStructure() );
  ss_profile->read_psipred_ss2(ss_file);

  std::string ss_as_string;
  for ( core::Size i = 1; i <= ss_profile->total_residue(); i++ ) {
    ss_as_string += ss_profile->secstruct(i);
  }

  std::cout << "ss " << ss_as_string << std::endl;

  std::string ss = ss_as_string;
  int pose_size = ipose.size() - 1;
  for (int i = 1; i < pose_size ; ++i) {
    ipose.set_secstruct(i, ss[i - 1]);
  }

  std::cout << ipose.secstruct() << std::endl;
  return ss_as_string;
}


inline void straigh_pose( core::pose::Pose& pose ) {
  for (int i = 1; i < static_cast<int>(pose.size()); ++i) {
    pose.set_phi(i, 0.0);
    pose.set_psi(i, 0.0);
  }

}

inline void read_pose(std::string pdb_path, core::pose::Pose& pose_ ) {
  using namespace protocols::simple_moves;
  using namespace core;

  std::cout << "pose from file" << std::endl;
  std::cout << "pdb_path " << pdb_path.c_str() << std::endl;
  //core::import_pose::pose_from_file(pose_, pdb_path.c_str(), false, core::import_pose::FileType::PDB_file);
  core::import_pose::centroid_pose_from_pdb(pose_, pdb_path.c_str(), false);
  std::string sequence = pose_.sequence();

  std::cout << "size seq " << pose_.size() << std::endl;
  // char z = 'Z';
  // int i = 1;
  // int size = pose_.size();

// if (pdb_path == "./input.pdb") {
//    pose_.delete_residue_range_slow(398, 400);
//   pose_.delete_residue_range_slow(795, 797);
//   pose_.delete_residue_slow(1193);
//   pose_.delete_residue_slow(1193);
//   pose_.delete_residue_slow(1192);
//   std::cout << "size seq " << pose_.size() << std::endl;
//  }
   // straigh_pose(ipose);
  std::cout << pose_.sequence() << std::endl;

  protocols::simple_moves::SwitchResidueTypeSetMover to_cen("centroid");
  to_cen.apply(pose_);
}

inline void read_density_map(std::string map_file, core::pose::Pose& ipose) {
  using namespace core::scoring::electron_density;
  ElectronDensity& densityMap = getDensityMap(map_file, true);
  std::cout << "start density align" << std::endl;
  protocols::electron_density::SetupForDensityScoringMoverOP dockindens( new protocols::electron_density::SetupForDensityScoringMover );
  dockindens->apply( ipose );
}


inline void read_frag_lib(std::string path, core::fragment::FragSetOP& frag_set_) {
  using namespace core::fragment;

  frag_set_ = FragmentIO().read_data( path.c_str() );
}


#endif // UTILS_HPP
