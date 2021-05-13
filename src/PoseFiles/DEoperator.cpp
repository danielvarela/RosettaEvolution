
#include <map>
#include <vector>
#include <string>
#include <stdlib.h>     /* exit, EXIT_FAILURE **/

#include "DEoperator.hpp"
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>

#include "../mpi_files/MasterRosettaCalculator.hpp"


DE_Operator::DE_Operator() {
  init_protein_archive();

  // score per stage
  score_per_stage["stage1"] = "score0";
  score_per_stage["stage2"] = "score1";
  score_per_stage["stage3"] = "score2";
  score_per_stage["stage4"] = "score3";

  // gmax_per_stage;
  gmax_per_stage["stage1"] = 30;
  gmax_per_stage["stage2"] = 10;
  gmax_per_stage["stage3"] = 10;
  gmax_per_stage["stage4"] = 10;

  // gmax_per_stage["stage1"] = 30;
  // gmax_per_stage["stage2"] = 10;
  // gmax_per_stage["stage3"] = 190;
  // gmax_per_stage["stage4"] = 100000;

  distances_map["error"] = dist_error;
  distances_map["rmsd"] = rmsd;
  distances_map["partial_rmsd"] = partial_rmsd;
  distances_map["rmsd_native_diff"] = rmsd_native_diff;
  distances_map["euclidean"] = euclidean;
  distances_map["euclidean_loop"] = euclidean_loop;
  distances_map["euclidean_diff_abs"] = euclidean_diff_abs;
  distances_map["euclidean_mario"] = euclidean_mario;
  distances_map["euclidean_partial_mario"] = euclidean_partial_mario;
  distances_map["mario_first_last"] = mario_first_last;

  protocol_name_map["error"] = protocol_error;
  protocol_name_map["Shared"] = Shared;
  protocol_name_map["HybridShared"] = HybridShared;
  protocol_name_map["HybridMover"] = HybridMover;
  protocol_name_map["CrowdingDE"] = CrowdingDE;
  protocol_name_map["MPICrowdingDE"] = MPICrowdingDE;
  protocol_name_map["MPISeedsDE"] = MPISeedsDE;
  protocol_name_map["MPIResetOldCrowdingDE"] = MPIResetOldCrowdingDE;
  protocol_name_map["HybridCrowdingDE"] = HybridCrowdingDE;
  protocol_name_map["ResetOld"] = ResetOldIndsHybrid;
  protocol_name_map["ResetOldCDE"] = ResetOldIndsCrowdingHybrid;

  init_popul_strategy_map["error"] = popul_error;
  init_popul_strategy_map["total_random"] = total_random;
  init_popul_strategy_map["random_pose_based"] = random_pose_based;
  init_popul_strategy_map["total_random_pose_based"] = total_random_pose_based;
  init_popul_strategy_map["init_popul_with_stage"] = init_popul_with_stage;

  fragment_insertion_strategy_map["error"] = frag_error;
  fragment_insertion_strategy_map["my_frag_insertion"] = my_frag_insertion;
  fragment_insertion_strategy_map["greedy_search"] = greedy_search;
  fragment_insertion_strategy_map["no_greedy_search"] = no_greedy_search;
  fragment_insertion_strategy_map["hybrid_mover"] = hybrid_mover;
  fragment_insertion_strategy_map["stage_rosetta_mover"] = stage_rosetta_mover;
  fragment_insertion_strategy_map["ILS_as_julia"] = ILS_as_julia;
}


void
DE_Operator::init_available_stages() {
  boost::char_separator<char> sep(",");
  std::string input_stages = app_options.get<std::string>("Protocol.stages");
  trim(input_stages);
  boost::tokenizer<boost::char_separator<char>> tokens(input_stages, sep);
  for (boost::tokenizer<boost::char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it) {
    vec_stages.push_back(*it);
  }
}

DE_Operator::DE_Operator(std::string prot_selected, boost::property_tree::ptree opt) : DE_Operator() {
  app_options = opt;
  fit_radius = app_options.get<double>("Extra.fitrad");
  init_files(prot_selection[prot_selected]);
  init_setup();
  init_two_stages_mover();
  init_available_stages();
}

DE_Operator::DE_Operator(std::string prot_selected) : DE_Operator() {
  init_files(prot_selection[prot_selected]);
  init_setup();
  init_default_score();
  init_popul = initialize_init_popul_strategy(std::string("total_random_pose_based"));
}


void
DE_Operator::init_protein_archive() {
  // Protinfo sap;
  // sap.name = "1sap";
  // sap.pdb_file = "./1sap_A.pdb";
  // sap.map_file = "./4054.map";
  // sap.frag_3 = "./frag3";
  // sap.frag_9 = "./frag9";
  // sap.ss_file = "./1sapA.ss2";

  Protinfo j9s;
  j9s.name = "3j9s";
  j9s.pdb_file = "./input.pdb";
  j9s.map_file = "./emd_6272.map";
  j9s.ss_file = "./sec.ss2";

  Protinfo c8c;
  c8c.name = "1c8cA";
  c8c.pdb_file = "./input_files/info_1c8cA/vf_1c8c.pdb";
  c8c.map_file = "./emd_6272.map";
  c8c.frag_3 = "./input_files/info_1c8cA/boinc_vf_aa1c8cA03_05.200_v1_3";
  c8c.frag_9 = "./input_files/info_1c8cA/boinc_vf_aa1c8cA09_05.200_v1_3";
  c8c.ss_file = "./input_files/info_1c8cA/vf_1c8cA.psipred_ss2";

  Protinfo c9o;
  c9o.name = "1c9oA";
  c9o.pdb_file = "./input_files/info_1c9oA/vf_1c9o.pdb";
  c9o.map_file = "./emd_6272.map";
  c9o.frag_3 = "./input_files/info_1c9oA/boinc_vf_aa1c9oA03_05.200_v1_3";
  c9o.frag_9 = "./input_files/info_1c9oA/boinc_vf_aa1c9oA09_05.200_v1_3";
  c9o.ss_file = "./input_files/info_1c9oA/vf_1c9oA.psipred_ss2";

  Protinfo ten;
  ten.name = "1ten";
  ten.pdb_file = "./input_files/info_1ten/vf_1ten.pdb";
  ten.map_file = "./emd_6272.map";
  ten.frag_3 = "./input_files/info_1ten/boinc_vf_aa1ten03_05.200_v1_3";
  ten.frag_9 = "./input_files/info_1ten/boinc_vf_aa1ten09_05.200_v1_3";
  ten.ss_file = "./input_files/info_1ten/vf_1ten.psipred_ss2";

  Protinfo bA;
  bA.name = "256bA";
  bA.pdb_file = "./input_files/info_256bA/vf_256bA.pdb";
  bA.map_file = "./emd_6272.map";
  bA.frag_3 = "./input_files/info_256bA/boinc_vf_aa256bA03_05.200_v1_3";
  bA.frag_9 = "./input_files/info_256bA/boinc_vf_aa256bA09_05.200_v1_3";
  bA.ss_file = "./input_files/info_256bA/vf_256bA.psipred_ss2";

  Protinfo elw;
  elw.name = "1elwA";
  elw.pdb_file = "./input_files/info_1elwA/vf_1elw.pdb";
  elw.map_file = "./emd_6272.map";
  elw.frag_3 = "./input_files/info_1elwA/boinc_vf_aa1elwA03_05.200_v1_3";
  elw.frag_9 = "./input_files/info_1elwA/boinc_vf_aa1elwA09_05.200_v1_3";
  elw.ss_file = "./input_files/info_1elwA/vf_1elwA.psipred_ss2";

  Protinfo opd;
  opd.name = "1opd";
  opd.pdb_file = "./input_files/info_1opd/vf_1opd.pdb";
  opd.map_file = "./emd_6272.map";
  opd.frag_3 = "./input_files/info_1opd/boinc_vf_aa1opd_03_05.200_v1_3";
  opd.frag_9 = "./input_files/info_1opd/boinc_vf_aa1opd_09_05.200_v1_3";
  opd.ss_file = "./input_files/info_1opd/vf_1opd_.psipred_ss2";

  Protinfo rnbA;
  rnbA.name = "1rnbA";
  rnbA.pdb_file = "./input_files/info_1rnbA/vf_1rnb.pdb";
  rnbA.map_file = "./emd_6272.map";
  rnbA.frag_3 = "./input_files/info_1rnbA/boinc_vf_aa1rnbA03_05.200_v1_3";
  rnbA.frag_9 = "./input_files/info_1rnbA/boinc_vf_aa1rnbA09_05.200_v1_3";
  rnbA.ss_file = "./input_files/info_1rnbA/vf_1rnbA.psipred_ss2";

  Protinfo fna;
  fna.name = "1fna";
  fna.pdb_file = "./input_files/info_1fna/vf_1fna.pdb";
  fna.map_file = "./emd_6272.map";
  fna.frag_3 = "./input_files/info_1fna/boinc_vf_aa1fna_03_05.200_v1_3";
  fna.frag_9 = "./input_files/info_1fna/boinc_vf_aa1fna_09_05.200_v1_3";
  fna.ss_file = "./input_files/info_1fna/vf_1fna_.psipred_ss2";

  Protinfo bgf;
  bgf.name = "1bgf";
  bgf.pdb_file = "./input_files/info_1bgf/vf_1bgf.pdb";
  bgf.map_file = "./emd_6272.map";
  bgf.frag_3 = "./input_files/info_1bgf/boinc_vf_aa1bgf_03_05.200_v1_3";
  bgf.frag_9 = "./input_files/info_1bgf/boinc_vf_aa1bgf_09_05.200_v1_3";
  bgf.ss_file = "./input_files/info_1bgf/vf_1bgf_.psipred_ss2";

  Protinfo who;
  who.name = "1who";
  who.pdb_file = "./input_files/info_1who/vf_1who.pdb";
  who.map_file = "./emd_6272.map";
  who.frag_3 = "./input_files/info_1who/boinc_vf_aa1who_03_05.200_v1_3";
  who.frag_9 = "./input_files/info_1who/boinc_vf_aa1who_09_05.200_v1_3";
  who.ss_file = "./input_files/info_1who/vf_1who_.psipred_ss2";

  Protinfo kpeA;
  kpeA.name = "1kpeA";
  kpeA.pdb_file = "./input_files/info_1kpeA/vf_1kpe.pdb";
  kpeA.map_file = "./emd_6272.map";
  kpeA.frag_3 = "./input_files/info_1kpeA/boinc_vf_aa1kpeA03_05.200_v1_3";
  kpeA.frag_9 = "./input_files/info_1kpeA/boinc_vf_aa1kpeA09_05.200_v1_3";
  kpeA.ss_file = "./input_files/info_1kpeA/vf_1kpeA.psipred_ss2";

  Protinfo dtdB;
  dtdB.name = "1dtdB";
  dtdB.pdb_file = "./input_files/info_1dtdB/1dtdB.pdb";
  dtdB.map_file = "./emd_6272.map";
  dtdB.frag_3 = "./input_files/info_1dtdB/aat000_03_05.200_v1_3";
  dtdB.frag_9 = "./input_files/info_1dtdB/aat000_09_05.200_v1_3";
  dtdB.ss_file = "./input_files/info_1dtdB/1dtdB.psipred_ss2";

  Protinfo sap;
  sap.name = "1sapA";
  sap.pdb_file = "./input_files/info_1sap/1sap.pdb";
  sap.map_file = "./emd_6272.map";
  sap.frag_3 = "./input_files/info_1sap/frag3";
  sap.frag_9 = "./input_files/info_1sap/frag9";
  sap.ss_file = "./input_files/info_1sap/1sapA.ss2";

  Protinfo wapA;
  wapA.name = "1wapA";
  wapA.pdb_file = "./input_files/info_1wapA/1wapA.pdb";
  wapA.map_file = "./emd_6272.map";
  wapA.frag_3 = "./input_files/info_1wapA/aat000_03_05.200_v1_3";
  wapA.frag_9 = "./input_files/info_1wapA/aat000_09_05.200_v1_3";
  wapA.ss_file = "./input_files/info_1wapA/1wapA.psipred_ss2";

  Protinfo ail;
  ail.name = "1ail";
  ail.pdb_file = "./input_files/info_1ail/1AIL.pdb";
  ail.map_file = "./emd_6272.map";
  ail.frag_3 = "./input_files/info_1ail/aa1ail03_05.200_v1_3";
  ail.frag_9 = "./input_files/info_1ail/aa1ail09_05.200_v1_3";
  ail.ss_file = "./input_files/info_1ail/1ail.psipred_ss2";






  Protinfo ail_nueva;
  ail_nueva.name = "1ail_nueva";
  ail_nueva.pdb_file = "./input_files/1ail_oldfrags_inputs/vf_1ail.pdb";
  ail_nueva.map_file = "./emd_6272.map";
  ail_nueva.frag_3 = "./input_files/1ail_oldfrags_inputs/boinc_vf_aa1ail_03_05.200_v1_3";
  ail_nueva.frag_9 = "./input_files/1ail_oldfrags_inputs/boinc_vf_aa1ail_09_05.200_v1_3";
  ail_nueva.ss_file = "./input_files/1ail_oldfrags_inputs/vf_1ail_.psipred_ss2";

  Protinfo hz6A;
  hz6A.name = "1hz6A";
  hz6A.pdb_file = "./input_files/1hz6A_oldfrags_inputs/vf_1hz6.pdb";
  hz6A.map_file = "./emd_6272.map";
  hz6A.frag_3 = "./input_files/1hz6A_oldfrags_inputs/boinc_vf_aa1hz6A03_05.200_v1_3";
  hz6A.frag_9 = "./input_files/1hz6A_oldfrags_inputs/boinc_vf_aa1hz6A09_05.200_v1_3";
  hz6A.ss_file = "./input_files/1hz6A_oldfrags_inputs/vf_1hz6A.psipred_ss2";

  Protinfo npsA;
  npsA.name = "1npsA";
  npsA.pdb_file = "./input_files/1npsA_oldfrags_inputs/vf_1npbs.pdb";
  npsA.map_file = "./emd_6272.map";
  npsA.frag_3 = "./input_files/1npsA_oldfrags_inputs/boinc_vf_aa1npsA03_05.200_v1_3";
  npsA.frag_9 = "./input_files/1npsA_oldfrags_inputs/boinc_vf_aa1npsA09_05.200_v1_3";
  npsA.ss_file = "./input_files/1npsA_oldfrags_inputs/vf_1npsA.psipred_ss2";

  Protinfo tig;
  tig.name = "1tig";
  tig.pdb_file = "./input_files/info_1tig/vf_1tig.pdb";
  tig.map_file = "./emd_6272.map";
  tig.frag_3 = "./input_files/info_1tig/boinc_vf_aa1tig_03_05.200_v1_3";
  tig.frag_9 = "./input_files/info_1tig/boinc_vf_aa1tig_09_05.200_v1_3";
  tig.ss_file = "./input_files/info_1tig/vf_1tig_.psipred_ss2";
  Protinfo tit;
  tit.name = "1tit";
  tit.pdb_file = "./input_files/info_1tit/vf_1tit.pdb";
  tit.map_file = "./emd_6272.map";
  tit.frag_3 = "./input_files/info_1tit/boinc_vf_aa1tit_03_05.200_v1_3";
  tit.frag_9 = "./input_files/info_1tit/boinc_vf_aa1tit_09_05.200_v1_3";
  tit.ss_file = "./input_files/info_1tit/vf_1tit_.psipred_ss2";

  Protinfo chf;
  chf.name = "2chf";
  chf.pdb_file = "./input_files/info_2chf/vf_2chf.pdb";
  chf.map_file = "./emd_6272.map";
  chf.frag_3 = "./input_files/info_2chf/boinc_vf_aa2chf_03_05.200_v1_3";
  chf.frag_9 = "./input_files/info_2chf/boinc_vf_aa2chf_09_05.200_v1_3";
  chf.ss_file = "./input_files/info_2chf/vf_2chf_.psipred_ss2";

  Protinfo acf;
  acf.name = "1acf";
  acf.pdb_file = "./input_files/info_1acf/vf_1acf.pdb";
  acf.map_file = "./emd_6272.map";
  acf.frag_3 = "./input_files/info_1acf/boinc_vf_aa1acf_03_05.200_v1_3";
  acf.frag_9 = "./input_files/info_1acf/boinc_vf_aa1acf_09_05.200_v1_3";
  acf.ss_file = "./input_files/info_1acf/vf_1acf_.psipred_ss2";
  Protinfo aiu;
  aiu.name = "1aiu";
  aiu.pdb_file = "./input_files/info_1aiu/vf_1aiu.pdb";
  aiu.map_file = "./emd_6272.map";
  aiu.frag_3 = "./input_files/info_1aiu/boinc_vf_aa1aiu_03_05.200_v1_3";
  aiu.frag_9 = "./input_files/info_1aiu/boinc_vf_aa1aiu_09_05.200_v1_3";
  aiu.ss_file = "./input_files/info_1aiu/vf_1aiu_.psipred_ss2";
  Protinfo a32;
  a32.name = "1a32";
  a32.pdb_file = "./input_files/info_1a32/vf_1a32.pdb";
  a32.map_file = "./emd_6272.map";
  a32.frag_3 = "./input_files/info_1a32/boinc_vf_aa1a32_03_05.200_v1_3";
  a32.frag_9 = "./input_files/info_1a32/boinc_vf_aa1a32_09_05.200_v1_3";
  a32.ss_file = "./input_files/info_1a32/vf_1a32_.psipred_ss2";
  Protinfo bk2;
  bk2.name = "1bk2";
  bk2.pdb_file = "./input_files/info_1bk2/vf_1bk2.pdb";
  bk2.map_file = "./emd_6272.map";
  bk2.frag_3 = "./input_files/info_1bk2/boinc_vf_aa1bk2_03_05.200_v1_3";
  bk2.frag_9 = "./input_files/info_1bk2/boinc_vf_aa1bk2_09_05.200_v1_3";
  bk2.ss_file = "./input_files/info_1bk2/vf_1bk2_.psipred_ss2";
  Protinfo a19;
  a19.name = "1a19A";
  a19.pdb_file = "./input_files/info_1a19/vf_1a19.pdb";
  a19.map_file = "./emd_6272.map";
  a19.frag_3 = "./input_files/info_1a19/boinc_vf_aa1a19A03_05.200_v1_3";
  a19.frag_9 = "./input_files/info_1a19/boinc_vf_aa1a19A09_05.200_v1_3";
  a19.ss_file = "./input_files/info_1a19/vf_1a19A.psipred_ss2";
  Protinfo b3a;
  b3a.name = "1b3aA";
  b3a.pdb_file = "./input_files/info_1b3a/vf_1b3a.pdb";
  b3a.map_file = "./emd_6272.map";
  b3a.frag_3 = "./input_files/info_1b3a/boinc_vf_aa1b3aA03_05.200_v1_3";
  b3a.frag_9 = "./input_files/info_1b3a/boinc_vf_aa1b3aA09_05.200_v1_3";
  b3a.ss_file = "./input_files/info_1b3a/vf_1b3aA.psipred_ss2";
  Protinfo bkr;
  bkr.name = "1bkrA";
  bkr.pdb_file = "./input_files/info_1bkr/vf_1bkr.pdb";
  bkr.map_file = "./emd_6272.map";
  bkr.frag_3 = "./input_files/info_1bkr/boinc_vf_aa1bkrA03_05.200_v1_3";
  bkr.frag_9 = "./input_files/info_1bkr/boinc_vf_aa1bkrA09_05.200_v1_3";
  bkr.ss_file = "./input_files/info_1bkr/vf_1bkrA.psipred_ss2";


  // 5 set of proteins
  Protinfo cg5;
  cg5.name = "1cg5B";
  cg5.pdb_file = "./input_files/info_1cg5B/vf_1cg5.pdb";
  cg5.map_file = "./emd_6272.map";
  cg5.frag_3 = "./input_files/info_1cg5B/boinc_vf_aa1cg5B03_05.200_v1_3";
  cg5.frag_9 = "./input_files/info_1cg5B/boinc_vf_aa1cg5B09_05.200_v1_3";
  cg5.ss_file = "./input_files/info_1cg5B/vf_1cg5B.psipred_ss2";
  Protinfo iib;
  iib.name = "1iibA";
  iib.pdb_file = "./input_files/info_1iibA/vf_1iib.pdb";
  iib.map_file = "./emd_6272.map";
  iib.frag_3 = "./input_files/info_1iibA/boinc_vf_aa1iibA03_05.200_v1_3";
  iib.frag_9 = "./input_files/info_1iibA/boinc_vf_aa1iibA09_05.200_v1_3";
  iib.ss_file = "./input_files/info_1iibA/vf_1iibA.psipred_ss2";


  Protinfo ctf;
  ctf.name = "1ctf";
  ctf.pdb_file = "./input_files/info_1ctf/vf_1ctf.pdb";
  ctf.map_file = "./emd_6272.map";
  ctf.frag_3 = "./input_files/info_1ctf/1ctf_.200.3mers";
  ctf.frag_9 = "./input_files/info_1ctf/1ctf_.200.9mers";
  ctf.ss_file = "./input_files/info_1ctf/1ctf_.psipred_ss2";

  Protinfo dhn;
  dhn.name = "1dhn";
  dhn.pdb_file = "./input_files/info_1dhn/vf_1dhn.pdb";
  dhn.map_file = "./emd_6272.map";
  dhn.frag_3 = "./input_files/info_1dhn/boinc_vf_aa1dhn_03_05.200_v1_3";
  dhn.frag_9 = "./input_files/info_1dhn/boinc_vf_aa1dhn_09_05.200_v1_3";
  dhn.ss_file = "./input_files/info_1dhn/vf_1dhn_.psipred_ss2";
  Protinfo lis;
  lis.name = "1lis";
  lis.pdb_file = "./input_files/info_1lis/vf_1lis.pdb";
  lis.map_file = "./emd_6272.map";
  lis.frag_3 = "./input_files/info_1lis/boinc_vf_aa1lis_03_05.200_v1_3";
  lis.frag_9 = "./input_files/info_1lis/boinc_vf_aa1lis_09_05.200_v1_3";
  lis.ss_file = "./input_files/info_1lis/vf_1lis_.psipred_ss2";

  Protinfo tul;
  tul.name = "1tul";
  tul.pdb_file = "./input_files/info_1tul/vf_1tul.pdb";
  tul.map_file = "./emd_6272.map";
  tul.frag_3 = "./input_files/info_1tul/boinc_vf_aa1tul_03_05.200_v1_3";
  tul.frag_9 = "./input_files/info_1tul/boinc_vf_aa1tul_09_05.200_v1_3";
  tul.ss_file = "./input_files/info_1tul/vf_1tul_.psipred_ss2";

  Protinfo vcc;
  vcc.name = "1vcc";
  vcc.pdb_file = "./input_files/info_1vcc/vf_1vcc.pdb";
  vcc.map_file = "./emd_6272.map";
  vcc.frag_3 = "./input_files/info_1vcc/boinc_vf_aa1vcc_03_05.200_v1_3";
  vcc.frag_9 = "./input_files/info_1vcc/boinc_vf_aa1vcc_09_05.200_v1_3";
  vcc.ss_file = "./input_files/info_1vcc/vf_1vcc_.psipred_ss2";

  Protinfo wit;
  wit.name = "1wit";
  wit.pdb_file = "./input_files/info_1wit/vf_1wit.pdb";
  wit.map_file = "./emd_6272.map";
  wit.frag_3 = "./input_files/info_1wit/boinc_vf_aa1wit_03_05.200_v1_3";
  wit.frag_9 = "./input_files/info_1wit/boinc_vf_aa1wit_09_05.200_v1_3";
  wit.ss_file = "./input_files/info_1wit/vf_1wit_.psipred_ss2";
  Protinfo vik;
  vik.name = "2vik";
  vik.pdb_file = "./input_files/info_2vik/vf_2vik.pdb";
  vik.map_file = "./emd_6272.map";
  vik.frag_3 = "./input_files/info_2vik/boinc_vf_aa2vik_03_05.200_v1_3";
  vik.frag_9 = "./input_files/info_2vik/boinc_vf_aa2vik_09_05.200_v1_3";
  vik.ss_file = "./input_files/info_2vik/vf_2vik_.psipred_ss2";

  Protinfo ci2;
  ci2.name = "2ci2I";
  ci2.pdb_file = "./input_files/info_2ci2I/vf_2ci2.pdb";
  ci2.map_file = "./emd_6272.map";
  ci2.frag_3 = "./input_files/info_2ci2I/boinc_vf_aa2ci2I03_05.200_v1_3";
  ci2.frag_9 = "./input_files/info_2ci2I/boinc_vf_aa2ci2I09_05.200_v1_3";
  ci2.ss_file = "./input_files/info_2ci2I/vf_2ci2I.psipred_ss2";


  // Prot5
  prot_selection["1cg5B"] = cg5;
  prot_selection["1ctf"] = ctf;
  prot_selection["1dhn"] = dhn;
  prot_selection["1iibA"] = iib;
  prot_selection["1lis"] = lis;
  prot_selection["1tul"] = tul;
  prot_selection["1vcc"] = vcc;
  prot_selection["1wit"] = wit;
  prot_selection["2vik"] = vik;
  prot_selection["2ci2I"] = ci2;

  prot_selection["1acf"] = acf;
  prot_selection["1bk2"] = bk2;
  prot_selection["1a19A"] = a19;
  prot_selection["1b3aA"] = b3a;
  prot_selection["1a32"] = a32;
  prot_selection["1aiu"] = aiu;
  prot_selection["1bkrA"] = bkr;


  prot_selection["1c8cA"] = c8c;
  prot_selection["1c9oA"] = c9o;
  prot_selection["1elwA"] = elw;
  prot_selection["1opd"] = opd;
  prot_selection["1rnbA"] = rnbA;
  prot_selection["1fna"] = fna;
  prot_selection["1who"] = who;
  prot_selection["1kpeA"] = kpeA;
  prot_selection["1chf"] = chf;
  prot_selection["1ten"] = ten;
  prot_selection["256bA"] = bA;

  // proteinas viejas
  prot_selection["1dtdB"] = dtdB;
  prot_selection["1sap"] = sap;
  prot_selection["1wapA"] = wapA;
  prot_selection["1ail"] = ail;
  prot_selection["1ail_nueva"] = ail_nueva;
  prot_selection["1hz6A"] = hz6A;

  prot_selection["1npsA"] = npsA;
  prot_selection["1tig"] = tig;
  prot_selection["1tit"] = tit;
  prot_selection["2chf"] = chf;

}

void
DE_Operator::init_setup() {
  pose_size = ipose.size();
  std::cout << "pose size " << pose_size << std::endl;
  pose_ = core::pose::PoseOP(new core::pose::Pose(ipose));
  best_pose_ = pose_->clone();
  native_pose_ = pose_->clone();
  scorefxn = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score3").c_str());

#if(USE_CRYO_EM)

  scorefxn->set_weight( core::scoring::elec_dens_fast, 30.0 );
  #endif

  std::cout << "start score: " << (*scorefxn)(*pose_) << std::endl;
  de_len = 10;
  start_res = 2;
  best_score = 1000000;
  // COMPLETE ABINITIO MOVER AT THE BOTTOM OF THE FILE
}

CalculateDistancePopulationPtr
DE_Operator::use_distances_strategy(std::string option) {
  return CalculateRmsdDistancePopulation(native_pose_, ffxn, ss, scorefxn, fit_radius).use_distances_strategy(option);
}


void
DE_Operator::init_files(Protinfo prot_info) {
  std::cout << "init protein " << prot_info.name << std::endl;
  read_pose(prot_info.pdb_file, ipose);
  ss = read_ss2(prot_info.ss_file, ipose);
  std::cout << "start read fragment lib " << std::endl;
  read_frag_lib(prot_info.frag_3, frag_set_ );
  read_frag_lib(prot_info.frag_9, frag_set_large );
  frag_opt.frag_set_ = frag_set_;
  frag_opt.frag_set_large = frag_set_large;


  core::pose::Pose my_pose;
  core::import_pose::centroid_pose_from_pdb(my_pose, "./out1_align.pdb", false);
  core::pose::PoseOP model_tmp = core::pose::PoseOP(new core::pose::Pose(my_pose));
  frag_opt.native_model = model_tmp;
  read_density_map("map_1wit.mrc", ipose);
}

boost::shared_ptr<MoverDE>
DE_Operator::prepare_stage(std::string stage_name) {
  /* INIT all funcions and classes neeeded by Differential Evolution
   */
  scorefxn = core::scoring::ScoreFunctionFactory::create_score_function(std::string(score_per_stage[stage_name]).c_str());
  frag_opt.scorefxn = scorefxn;
  frag_opt.stage_name = stage_name;
  core::pose::Pose my_pose;
  core::import_pose::centroid_pose_from_pdb(my_pose, "./out1_align.pdb", false);
  core::pose::PoseOP model_tmp = core::pose::PoseOP(new core::pose::Pose(my_pose));
  frag_opt.native_model = model_tmp;

  bool fragment_at_trials = false;
  boost::property_tree::ptree::const_assoc_iterator it = app_options.find("Fragments.strategy_at_trials");
  if(app_options.not_found() == it) {
    fragment_at_trials = true;
    frag_mover = initialize_fragment_insertion_strategy(app_options.get<std::string>("Fragments.strategy_at_trials"));
  } else {
    fragment_at_trials = false;
    frag_mover = initialize_fragment_insertion_strategy("my_frag_insertion");
  }

  if (fragment_at_trials) {
    ffxn = boost::shared_ptr<FitFunction>( new PoseFragmentFunction(pose_, scorefxn, ss, frag_mover));
#if(USE_CRYO_EM)
    ffxn = boost::shared_ptr<FitFunction>( new PoseDensityFragmentFunction(pose_, scorefxn, ss, frag_mover));
#endif

  } else {
    //    ffxn = boost::shared_ptr<FitFunction>( new PoseToyFunction(pose_, scorefxn, ss, frag_mover));
    ffxn = boost::shared_ptr<FitFunction>( new PoseScoreFunction(pose_, scorefxn, ss, frag_mover));
#if(USE_CRYO_EM)

    ffxn = boost::shared_ptr<FitFunction>( new PoseDensityFunction(pose_, scorefxn, ss, frag_mover));
  #endif


  }


   if ( (app_options.get<std::string>("Protocol.name") =="HybridShared") ||
	(app_options.get<std::string>("Protocol.name") =="MPIResetOldCrowdingDE") ||
	(app_options.get<std::string>("Protocol.name") =="HybridMover") ||
	(app_options.get<std::string>("Protocol.name") =="MPICrowdingDE") ||
	(app_options.get<std::string>("Protocol.name") =="MPISeedsDE") ||
	(app_options.get<std::string>("Protocol.name") =="ResetOld") ||
	(app_options.get<std::string>("Protocol.name") =="ResetOldCDE") ||
	(app_options.get<std::string>("Protocol.name") =="HybridCrowdingDE")) {
     bool fragment_at_popul = false;
     boost::property_tree::ptree::const_assoc_iterator it = app_options.find("Fragments.strategy_at_population");
     if(app_options.not_found() == it) {
       fragment_at_popul = true;
       initialize_local_search_to_apply_at_population(app_options.get<std::string>("Fragments.strategy_at_population"));
     } else {
       //default option
       initialize_local_search_to_apply_at_population(std::string("stage_rosetta_mover"));
     }
   }


  ffxn->name_ = std::string("Score Function for stage " + stage_name);
  print_best = PrintBestIndividualPtr( new PrintBestIndividual(pose_, ffxn, ss, scorefxn));



  //calculate_distances_popul = use_distances_strategy("rmsd");
  calculate_distances_popul = use_distances_strategy( app_options.get<std::string>("Protocol.distance_strategy") );
  update_current_population_to_stage_score();

  boost::shared_ptr<MoverDE> de = init_differential_evolution_protocol();
  de->Gmax = gmax_per_stage[stage_name];
  de->use_print_class = true;
  de->print_best_ind = print_best;
  de->calculate_distances_popul = calculate_distances_popul;
  CalculateDistancePopulationPtr calculate_partial_rmsd = use_distances_strategy( "partial_rmsd" );
  de->calculator_partial_rmsd = boost::dynamic_pointer_cast<PartialRMSDcalculator>(calculate_partial_rmsd);

  //  de->popul = current_population;

  return de;
}

void
DE_Operator::initialize_local_search_to_apply_at_population(std::string input_option) {
  local_search_frag_mover = initialize_fragment_insertion_strategy(input_option);
  local_search = boost::shared_ptr<LocalSearchIndividualMover>(new LocalSearchIndividualMover(pose_, scorefxn, ss, local_search_frag_mover));
}

boost::shared_ptr<FragInsertionMover>
DE_Operator::initialize_fragment_insertion_strategy(std::string input_option ) {
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
  case hybrid_mover: {
     strategy_return = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::hybrid_mover, frag_opt);
    break;
  }
  case ILS_as_julia: {
    strategy_return = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::ILS_as_julia, frag_opt);
    break;
  }
  default:
    std::cout << "no insertion strategy option" << std::endl;
    exit(1);
    break;
  }
  return strategy_return;


}

boost::shared_ptr<InitPopulation>
DE_Operator::initialize_init_popul_strategy(std::string init_popul_option) {
  boost::shared_ptr<InitPopulation> init_popul_in_options_file;

  switch (init_popul_strategy_map[init_popul_option]) {
  case total_random: {
    init_popul_in_options_file = boost::shared_ptr<InitPopulation>(new InitPopulation(ffxn));
    break;
  }
  case total_random_pose_based: {
    init_popul_in_options_file = boost::shared_ptr<InitPopulation>(new PoseTotalRandomInitPopulation(ffxn, pose_, ss));
    break;
  }
  case random_pose_based: {
    init_popul_in_options_file = boost::shared_ptr<InitPopulation>(new PosePopulation(ffxn, pose_, ss));
    break;
  }
  case init_popul_with_stage: {
    init_popul_in_options_file = boost::shared_ptr<InitPopulation>(new TwoStagesInitPopulation(ffxn, pose_, two_stages_mover, ss));
    break;
  }
  default:
    std::cout << "problem selecting the init population strategy" << std::endl;
    exit(1);
    break;
  }
  return init_popul_in_options_file;
}

void
DE_Operator::init_default_score() {
  core::scoring::ScoreFunctionOP score3 = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score3").c_str());
  frag_opt.scorefxn = score3;
  frag_opt.stage_name = "stage4";
  frag_mover = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::greedy_search, frag_opt);
  ffxn = boost::shared_ptr<FitFunction>( new PoseScoreFunction(pose_, score3, ss, frag_mover)); 
#if(USE_CRYO_EM)

    ffxn = boost::shared_ptr<FitFunction>( new PoseDensityFunction(pose_, scorefxn, ss, frag_mover));
  #endif


}

void
DE_Operator::init_two_stages_mover() {
  core::scoring::ScoreFunctionOP score1 = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score1").c_str());
  frag_opt.scorefxn = score1;
  frag_opt.stage_name = "stage1";
  frag_mover = FragInsertionStrategy::get(FragInsertionStrategy::FragMoverTypes::greedy_search, frag_opt);
  ffxn = boost::shared_ptr<FitFunction>( new PoseScoreFunction(pose_, score1, ss, frag_mover));
#if(USE_CRYO_EM)
    ffxn = boost::shared_ptr<FitFunction>( new PoseDensityFunction(pose_, scorefxn, ss, frag_mover));
  #endif


  two_stages_mover = boost::shared_ptr<InitStagesMover>( new InitStagesMover(score1, frag_set_, frag_set_large));
    init_popul = initialize_init_popul_strategy(app_options.get<std::string>("Protocol.init_strategy") );
}

void
DE_Operator::update_current_population_to_stage_score( ) {
  if (current_population.size() == 0) {
    //    current_population = de->popul;
    std::cout << "error, the current population was not initialize" << std::endl;
    exit(1);
  } else {
    CalculateDistancePopulation::DistancesResult result =
      calculate_distances_popul->run(current_population);
    std::vector<double> current_rmsd = result.rmsd_to_native;
    // rescore population
    for (int i = 0; i < current_population.size(); i++) {
      ffxn->score(current_population[i]);
    }
  }
}

boost::shared_ptr<MoverDE>
DE_Operator::init_differential_evolution_protocol() {
  //     boost::shared_ptr<MoverDE> de = boost::shared_ptr<MoverDE>(new SharedMoverDE(ffxn, init_popul));
  ConfigurationDE conf = ConfigurationDE(app_options.get<int>("DE.NP"), app_options.get<double>("DE.CR"), app_options.get<double>("DE.F"), app_options.get<double>("Extra.fitrad") , app_options.get<std::string>("Protocol.prot"));

  std::string protocol_name = app_options.get<std::string>("Protocol.name");
  std::string crossover_strategy = app_options.get<std::string>("Protocol.crossover_strategy");
  std::string select_parents_strategy = app_options.get<std::string>("Protocol.select_parents_strategy");

  boost::shared_ptr<MoverDE> de ;
  switch (protocol_name_map[protocol_name]) {
  case Shared: {
    de = boost::shared_ptr<MoverDE>(new SharedMoverDE(app_options, ffxn, current_population));
    break;
  }
  case HybridMover: {
    de = boost::shared_ptr<MoverDE>(new HybridMoverDE(app_options, ffxn, current_population, local_search));
    break;
  }
  case HybridShared: {
    de = boost::shared_ptr<MoverDE>(new SharedHybridMoverDE(app_options, ffxn, current_population, local_search));
    break;
  }
  case CrowdingDE: {
    de = boost::shared_ptr<MoverDE>(new CrowdingMoverDE(app_options, ffxn, current_population));
    break;
  }
#if(MPI_ENABLED)
  case MPIResetOldCrowdingDE: {
    de = boost::shared_ptr<MoverDE>(new MPIResetOldCrowdingHybridDE(app_options, ffxn, current_population, local_search, boost::shared_ptr<LocalSearchIndividualMover>(new LocalSearchIndividualMover(pose_, scorefxn, ss, initialize_fragment_insertion_strategy("ILS_as_julia")  ))  ));
    break;
  }
  case MPISeedsDE: {
    de = boost::shared_ptr<MoverDE>(new MPISeedsMoverDE(app_options, ffxn, current_population, local_search));
    break;
  }
  case MPICrowdingDE: {
    de = boost::shared_ptr<MoverDE>(new MPICrowdingMoverDE(app_options, ffxn, current_population, local_search));
    break;
  }
#endif
  case ResetOldIndsHybrid: {
    de = boost::shared_ptr<MoverDE>(new ResetOldIndsHybridDE(app_options, ffxn, current_population, local_search, boost::shared_ptr<LocalSearchIndividualMover>(new LocalSearchIndividualMover(pose_, scorefxn, ss, initialize_fragment_insertion_strategy("ILS_as_julia")  ))   ));
    break;
  }
  case ResetOldIndsCrowdingHybrid: {
    de = boost::shared_ptr<MoverDE>(new ResetOldIndsCrowdingHybridDE(app_options, ffxn, current_population, local_search, boost::shared_ptr<LocalSearchIndividualMover>(new LocalSearchIndividualMover(pose_, scorefxn, ss, initialize_fragment_insertion_strategy("ILS_as_julia")  ))   ));
    break;
  }
  case HybridCrowdingDE: {
    de = boost::shared_ptr<MoverDE>(new CrowdingHybridMoverDE(app_options, ffxn, current_population, local_search));
    break;
  }
default:
  std::cout << "error with Protocol Name" << std::endl;
  exit(1);
    break;
  }
  de->init_popul_ = initialize_init_popul_strategy("total_random_pose_based");
  de->select_parents_strategy = select_parents_strategy;
  de->crossover_strategy = crossover_strategy;

  return de;
}

void
DE_Operator::print_best_pose() {
  double obtained_score = (*scorefxn)(*best_pose_);
  std::cout << "obtained score: " << obtained_score << std::endl;
  for (int i = 1; i < static_cast<int>(best_pose_->size()); ++i) {
    std::cout << best_pose_->phi(i) << " " << best_pose_->psi(i) << " ";
  }
  std::cout << std::endl;
  best_pose_->dump_pdb("./result_relax.pdb");
}

void
DE_Operator::apply_differential_evolution_with_stage(std::string stage_name) {
  boost::shared_ptr<MoverDE> de;
  de = prepare_stage(stage_name);
  de->apply();
  current_population = de->popul;
  best_pose_ = pose_->clone();
  ffxn->fill_pose(best_pose_, de->best_ind(), ss);
  std::cout << "BEST POSE FOUND AT STAGE " << stage_name << " " << de->best_ind().score << std::endl;
}

core::pose::PoseOP
DE_Operator::get_native_pose() {
  return native_pose_;
}

void
DE_Operator::initialize_population_stage1() {
  // current_population is modified using the stage1 of rosetta
  std::cout << "=== init the current population ================ " << std::endl;
  init_popul->apply(current_population, app_options.get<int>("DE.NP"), ffxn->D() );
  std::cout << "=== finish initialization current population === " << std::endl;
}

void
DE_Operator::run_complete_abinitio() {
#if(MPI_ENABLED)


    core::scoring::ScoreFunctionOP score_4 = core::scoring::ScoreFunctionFactory::create_score_function(std::string("score3").c_str());
    boost::shared_ptr<CompleteAbinitioMover> abinitio = boost::shared_ptr<CompleteAbinitioMover>(new CompleteAbinitioMover(score_4, frag_set_, frag_set_large));
    prepare_stage("stage4");
    ffxn = boost::shared_ptr<FitFunction>( new PoseFragmentFunction(pose_, scorefxn, ss, abinitio));
#if(USE_CRYO_EM)
    ffxn = boost::shared_ptr<FitFunction>( new PoseDensityFragmentFunction(pose_, scorefxn, ss, frag_mover));
#endif
    ffxn->set_name("complete");
    boost::shared_ptr<FitFunction> simple_ffxn = boost::shared_ptr<FitFunction>( new PoseScoreFunction(pose_, scorefxn, ss, frag_mover));

#if(USE_CRYO_EM)
    simple_ffxn = boost::shared_ptr<FitFunction>( new PoseDensityFunction(pose_, scorefxn, ss, frag_mover));
#endif



  boost::shared_ptr<MasterRosettaCalculator> mpi_calculator;
  mpi_calculator = boost::shared_ptr<MasterRosettaCalculator>(new MasterRosettaCalculator(ffxn));
  current_population = mpi_calculator->run(current_population);
    std::cout << "finish rosetta" << std::endl;
    print_final_population(current_population);
#endif
}

void DE_Operator::print_final_population(std::vector<Individual> popul ) {
  CalculateDistancePopulation::DistancesResult result =
    calculate_distances_popul->run(popul);
  std::vector<double> current_rmsd = result.rmsd_to_native;
  std::cout << "[POP] ";
  double acc = 0;
  double best = 100000;
  for (int i = 0; i < app_options.get<int>("DE.NP"); i++) {
    double score =  (-1 * SCORE_ERROR_FIXED) +  popul[i].score;
    acc += score;
    if (score < best) {
      best = score;
    }
  }
  std::cout << "[GEN] 99999 " << best << " " << acc/popul.size() << std::endl;
   for (int i = 0; i < app_options.get<int>("DE.NP"); i++) {
    double score = popul[i].score;
    std::cout << (-1 * SCORE_ERROR_FIXED) +  score << " , ";
  }

  std::cout << std::endl;
  std::cout << "[RMSD_NATIVE] ";
  for (int i = 0; i < app_options.get<int>("DE.NP"); i++) {
    double rmsd = current_rmsd[i];
    std::cout << rmsd << " , ";
  }
  std::cout << std::endl;
}

void
DE_Operator::run() {
  std::string stage_name;
  std::cout << "run " << std::endl;

  std::string protocol_name = app_options.get<std::string>("Protocol.name");



  if (false) {
    initialize_population_stage1();
    run_complete_abinitio();
  } else {
    initialize_population_stage1();
    PrintRMSDvsMarioAnalisys printer_mario_pairs(  use_distances_strategy("rmsd"),  use_distances_strategy("euclidean_partial_mario"));
    std::cout << "================================================" << std::endl;
    printer_mario_pairs.apply(current_population);
    std::cout << "================================================" << std::endl;

    if ( std::find(vec_stages.begin(), vec_stages.end(), "stage2") != vec_stages.end()) {
      stage_name = "stage2";
      std::cout << "================================================" << std::endl;
      std::cout << "START STAGE " << stage_name << std::endl;
      apply_differential_evolution_with_stage(stage_name);
    }


    if ( std::find(vec_stages.begin(), vec_stages.end(), "stage3") != vec_stages.end()) {
      stage_name = "stage3";
      std::cout << "================================================" << std::endl;
      std::cout << "START STAGE " << stage_name << std::endl;
      apply_differential_evolution_with_stage(stage_name);
    }

    if ( std::find(vec_stages.begin(), vec_stages.end(), "stage4") != vec_stages.end()) {
      stage_name = "stage4";
      std::cout << "================================================" << std::endl;
      std::cout << "START STAGE " << stage_name << std::endl;
      apply_differential_evolution_with_stage(stage_name);
    }




    std::cout << "================================================" << std::endl;
    printer_mario_pairs.apply(current_population);
    std::cout << "================================================" << std::endl;
    //print_best_pose();
    print_final_population(current_population);
  }
}
