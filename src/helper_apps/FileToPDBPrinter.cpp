
#include "../inc/FileToPDBPrinter.hpp"

void
FileToPDBPrinter::population_to_file_json(boost::property_tree::ptree pt, std::string score_name, int gen, const std::vector<Individual>& popul) {
  int NP = pt.get<int>("DE.NP");
  int D = popul[0].vars.size();
  boost::char_separator<char> sep(" ");
  tokenizer tokens(score_name, sep);
  std::string sc_name;
  tokenizer::iterator it = tokens.begin();
  for (; it != tokens.end(); ++it) {
    sc_name = *it;
  }

  double fit_rad = pt.get<double>("Extra.fitrad");
  std::string prot = pt.get<std::string>("Protocol.prot");
  std::string output_folder = pt.get<std::string>("PrintPopul.output");

  std::ofstream ofs ( string("./"+output_folder+"/population_"+sc_name+"_"+ to_string(gen) +".log") , std::ofstream::out);
  //ofs << "[CONF] " << "NP " << NP << " D " << D << " score " << sc_name << " fit_rad "<<  fit_rad << " prot " << prot << " gen " << gen << std::endl;
  ofs << "{" << std::endl;
  ofs << " conf : { " << "NP : " << NP << " ,  D : " << D << " ,  score : " << sc_name << " ,  fit_rad : "<<  fit_rad << " , prot : " << prot << " ,  gen : " << gen << "} , " << std::endl;

  ofs << " individuals : [";
  for (size_t i = 0; i < popul.size(); i++) {
    ofs << " [ ";
    std::vector<double> vars = popul[i].vars;
    for (size_t j = 0; j < vars.size(); j++) {
      if (j == (vars.size() - 1)) {
	ofs << vars[j] << " ] ";
      } else {
	ofs << vars[j] << " , ";
      }
    }
    if (i == (popul.size() - 1)) {
      ofs << std::endl << " ] " << std::endl;
    } else {
      ofs << " , " << std::endl;
    }
  }
  ofs << "}" << std::endl;
  ofs.close();

}
