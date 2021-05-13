
#include "ProteinInfoArchive.hpp"

void
ProteinInfoArchive::build() {
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
  npsA.pdb_file = "./input_files/info_1npsA/vf_1nps.pdb";
  npsA.map_file = "./emd_6272.map";
  npsA.frag_3 = "./input_files/info_1npsA/boinc_vf_aa1npsA03_05.200_v1_3";
  npsA.frag_9 = "./input_files/info_1npsA/boinc_vf_aa1npsA09_05.200_v1_3";
  npsA.ss_file = "./input_files/info_1npsA/vf_1npsA.psipred_ss2";

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
  bkr.pdb_file = "./input_files/info_1bkrA/vf_1bkr.pdb";
  bkr.map_file = "./emd_6272.map";
  bkr.frag_3 = "./input_files/info_1bkrA/boinc_vf_aa1bkrA03_05.200_v1_3";
  bkr.frag_9 = "./input_files/info_1bkrA/boinc_vf_aa1bkrA09_05.200_v1_3";
  bkr.ss_file = "./input_files/info_1bkrA/vf_1bkrA.psipred_ss2";


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

  Protinfo dimaio;
  dimaio.name = "dimaio";
  dimaio.pdb_file = "./input_files/info_dimaio/starting_dimaio.pdb";
  dimaio.map_file = "./input_files/info_dimaio/monomer_nolig.mrc";
  dimaio.frag_3 = "./input_files/info_dimaio/5mers";
  dimaio.frag_9 = "./input_files/info_dimaio/5mers";
  dimaio.ss_file = "./input_files/info_dimaio/dimaio.psipred_ss2";
  Protinfo eyv;
  eyv.name = "1eyvA";
  eyv.pdb_file = "./input_files/info_1eyvA/vf_1eyv.pdb";
  eyv.map_file = "./emd_6272.map";
  eyv.frag_3 = "./input_files/info_1eyvA/boinc_vf_aa1eyvA03_05.200_v1_3";
  eyv.frag_9 = "./input_files/info_1eyvA/boinc_vf_aa1eyvA09_05.200_v1_3";
  eyv.ss_file = "./input_files/info_1eyvA/vf_1eyvA.psipred_ss2";

  Protinfo gvp;
  gvp.name = "1gvp";
  gvp.pdb_file = "./input_files/info_1gvp/vf_1gvp.pdb";
  gvp.map_file = "./emd_6272.map";
  gvp.frag_3 = "./input_files/info_1gvp/boinc_vf_aa1gvp_03_05.200_v1_3";
  gvp.frag_9 = "./input_files/info_1gvp/boinc_vf_aa1gvp_09_05.200_v1_3";
  gvp.ss_file = "./input_files/info_1gvp/vf_1gvp_.psipred_ss2";


  Protinfo nc8c;
  nc8c.name = "n1c8cA";
  nc8c.frag_3 =  "./input_files/1c8cA_newfrags_inputs/1c8cA.200.3mers";
  nc8c.frag_9 = "./input_files/1c8cA_newfrags_inputs/1c8cA.200.9mers";
  nc8c.ss_file =  "./input_files/1c8cA_newfrags_inputs/1c8cA.psipred_ss2";
  nc8c.pdb_file =  "./input_files/1c8cA_newfrags_inputs/vf_1c8c.pdb";
  nc8c.map_file =  "./input_files/1c8cA_newfrags_inputs/vf_1c8c.pdb";

  Protinfo nelw;

  nelw.name = "n1elwA";
  nelw.frag_3 =   "./input_files/1elwA_newfrags_inputs/1elwA.200.3mers";
  nelw.frag_9 =   "./input_files/1elwA_newfrags_inputs/1elwA.200.9mers";
  nelw.ss_file = "./input_files/1elwA_newfrags_inputs/1elwA.psipred_ss2";
  nelw.pdb_file =   "./input_files/1elwA_newfrags_inputs/vf_1elw.pdb";

  Protinfo nhz6;
  nhz6.name = "n1hz6A";
  nhz6.frag_3 =   "./input_files/1hz6A_newfrags_inputs/1hz6A.200.3mers";
  nhz6.frag_9 =   "./input_files/1hz6A_newfrags_inputs/1hz6A.200.9mers";
  nhz6.ss_file =   "./input_files/1hz6A_newfrags_inputs/1hz6A.psipred_ss2";
  nhz6.pdb_file =     "./input_files/1hz6A_newfrags_inputs/vf_1hz6.pdb";

  Protinfo nkpe;
  nkpe.name = "n1kpeA";
  nkpe.frag_3 =   "./input_files/1kpeA_newfrags_inputs/1kpeA.200.3mers";
  nkpe.frag_9 =   "./input_files/1kpeA_newfrags_inputs/1kpeA.200.9mers";
  nkpe.ss_file =     "./input_files/1kpeA_newfrags_inputs/1kpeA.psipred_ss2";
  nkpe.pdb_file =     "./input_files/1kpeA_newfrags_inputs/vf_1kpe.pdb";

  Protinfo nopd;
  nopd.name = "n1nopd";
  nopd.frag_3 =   "./input_files/1opd_newfrags_inputs/1opd_.200.3mers";
  nopd.frag_9 =   "./input_files/1opd_newfrags_inputs/1opd_.200.9mers";
  nopd.ss_file =     "./input_files/1opd_newfrags_inputs/1opd_.psipred_ss2";
  nopd.pdb_file =     "./input_files/1opd_newfrags_inputs/vf_1opd.pdb";

  Protinfo nwit;
  nwit.name = "n1wit";
  nwit.frag_3 =   "./input_files/1wit_newfrags_inputs/1wit_.200.3mers";
  nwit.frag_9 =   "./input_files/1wit_newfrags_inputs/1wit_.200.9mers";
  nwit.ss_file = "./input_files/1wit_newfrags_inputs/1wit_.psipred_ss2";
  nwit.pdb_file = "./input_files/1wit_newfrags_inputs/vf_1wit.pdb";
  
  prot_selection["n1c8cA"] = nc8c;
  prot_selection["n1elwA"] = nelw;
  prot_selection["n1hz6A"] = nhz6;
  prot_selection["n1kpeA"] = nkpe;
  prot_selection["n1opd"] = nopd;
  prot_selection["n1wit"] = nwit;
  
  // Prot5
  prot_selection["1cg5B"] = cg5;
  prot_selection["1bgf"] = bgf;
  prot_selection["1gvp"] = gvp;
  prot_selection["1eyvA"] = eyv;
  prot_selection["dimaio"] = dimaio;
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
