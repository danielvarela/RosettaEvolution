
#ifndef DISTANCESPROFILECALCULATOR_H
#define DISTANCESPROFILECALCULATOR_H

#include "DE_types.hpp"
#include "FitFunction.hpp"
#include "FragInsertionMover.hpp"
#include "PoseFunction.hpp"

typedef std::map<double, int> DistancesProfileMap;

class DistanceProfileRes
{
public:
  int index_;
  int index_2_;
  DistancesProfileMap map_;
  std::vector<double> bin_values;
  DistanceProfileRes(int ind1, int ind2) ;

  double scale( double old_value) ;

  void add_new_calculation(const Individual& ind1,const Individual& ind2) ;

  double get_peak() ;
};


class DistancesProfileCalculator
{
public:
  DistancesProfileCalculator() {}

  void apply(const std::vector<Individual>& popul);
};


#endif /* DISTANCESPROFILECALCULATOR_H */
