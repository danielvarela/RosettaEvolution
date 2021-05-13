
#include "FitFunction.hpp"
#include <vector>
#include <map>
#include <string>

double
FitFunction::scale( double old_value) {
  double new_min = lim() * -1;
  double new_max = lim();

  return ( (old_value - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min;
}

Individual
FitFunction::convert(const Individual& ind) {

  Individual ind_converted = ind;
  for (int i =0 ; i < D(); ++i) {
    ind_converted.vars[i] = scale(ind.vars[i]);
  }

  return ind_converted;
}
