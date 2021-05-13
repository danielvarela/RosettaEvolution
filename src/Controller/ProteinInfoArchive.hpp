#ifndef PROTEININFOARCHIVE_H
#define PROTEININFOARCHIVE_H


#include <map>
#include <vector>
#include <string>

#include "../Extra/Utils.hpp"

class ProteinInfoArchive
{
public:
  std::map<std::string, Protinfo> prot_selection;
  ProteinInfoArchive() {
    build();
  }

  void build();

  std::map<std::string, Protinfo> get_map() {
    return prot_selection;
  }
};


#endif /* PROTEININFOARCHIVE_H */
