#include "FragmentInserterBuilder.hpp"

FragInsertionStrategy::FragOptions
FragmentInserterBuilder::init_frag_files(Protinfo prot_info, boost::property_tree::ptree app_options ) {
  FragInsertionStrategy::FragOptions frag_opt;
  core::pose::PoseOP pose_, best_pose_, native_pose_;
  core::pose::Pose ipose;
  core::fragment::FragSetOP frag_set_, frag_set_large;
  read_pose(prot_info.pdb_file, ipose);

  
#if(USE_CRYO_EM)
  init_density_map(ipose, prot_info.name, app_options.get<int>("Protocol.map_res"));

#endif

  // int size_pose = ipose.total_residue();
  // if ( size_pose < ipose.total_residue()) {
  //   core::pose::Pose recortado = core::pose::Pose( ipose ,1, size_pose);
  //   ipose = recortado;
  // }
  // if (!ipose.is_centroid()) {
  //     protocols::moves::MoverOP tocen( new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID ) );
  //     tocen->apply(ipose);
  //   }

  core::pose::Pose my_pose;
  read_pose(prot_info.pdb_file, my_pose);

#if(USE_CRYO_EM)
  init_density_map(my_pose, prot_info.name, app_options.get<int>("Protocol.map_res"));
#endif

  std::string ss = read_ss2(prot_info.ss_file, ipose);
  //std::cout << "start read fragment lib " << std::endl;
  frag_opt.ss = ss;
  read_frag_lib(prot_info.frag_3, frag_set_ );
  read_frag_lib(prot_info.frag_9, frag_set_large );
  frag_opt.frag_set_ = frag_set_;
  frag_opt.frag_set_large = frag_set_large;
  // core::pose::Pose my_pose;
  // core::import_pose::centroid_pose_from_pdb(my_pose, "./input_files/out1_align.pdb", false);


  core::pose::PoseOP model_tmp = core::pose::PoseOP(new core::pose::Pose(my_pose));
  core::pose::PoseOP model_native = core::pose::PoseOP(new core::pose::Pose(ipose));
  frag_opt.native_model = model_native;
  frag_opt.consen_model = model_tmp;
  return frag_opt;
}

void
FragmentInserterBuilder::init_density_map(    core::pose::Pose& ipose, std::string prot_name, int resolution ) {
  switch (resolution) {
  case 3: {
    if (prot_name =="dimaio") {
      read_density_map("monomer_nolig.mrc", ipose);
    }  
  }
  case 5: {
    if (prot_name =="dimaio") {
      read_density_map("monomer_nolig.mrc", ipose);
    }  
    if (prot_name =="1wit") {
      read_density_map("map_1wit.mrc", ipose);
    }
    if (prot_name =="1c9oA") {
      read_density_map("map_1c9oA.mrc", ipose);
    }

    if (prot_name =="256bA") {
      read_density_map("map_256bA.mrc", ipose);
    }
    break;
  }
  case 2: {
    if (prot_name =="dimaio") {
      read_density_map("monomer_nolig.mrc", ipose);
    }  
    if (prot_name =="1wit") {
      read_density_map("map_1wit_2A.mrc", ipose);
    }
    if (prot_name =="1c9oA") {
      read_density_map("map_1c9oA_2A.mrc", ipose);
    }

    if (prot_name =="256bA") {
      read_density_map("map_256bA_2A.mrc", ipose);
    }
    break;
  }
  case 10: {
    if (prot_name =="dimaio") {
      read_density_map("monomer_nolig.mrc", ipose);
    } 
    if (prot_name =="1wit") {
      read_density_map("map_1wit_10.mrc", ipose);
    }
    if (prot_name =="1c9oA") {
      read_density_map("map_1c9oA_10A.mrc", ipose);
    }

    if (prot_name =="256bA") {
      read_density_map("map_256bA_10A.mrc", ipose);
    }
    break;
  }
default:
    break;
  }
}


boost::shared_ptr<FragInsertionMover>
FragmentInserterBuilder::get(std::string input_option, FragInsertionStrategy::FragOptions frag_opt) {
  boost::shared_ptr<FragInsertionMover> strategy_return;
  switch (fragment_insertion_strategy_map[input_option]) {
  case my_frag_insertion: {
    strategy_return = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::my_frag_insertion, frag_opt);
    break;
  }
  case greedy_search: {
    strategy_return = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::greedy_search, frag_opt);
    break;
  }
  case no_greedy_search: {
    strategy_return = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::no_greedy_search, frag_opt);
    break;
  }
  case stage_rosetta_mover: {
    strategy_return = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::stage_rosetta_mover, frag_opt);
    break;
  }
  default:
    std::cout << "no insertion strategy option" << std::endl;
    exit(1);
    break;
  }
  return strategy_return;
}
