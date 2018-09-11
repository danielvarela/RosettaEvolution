
#include "DistancesProfileCalculator.hpp"
#include <map>
#include <vector>


DistanceProfileRes::DistanceProfileRes(int ind1, int ind2) {
  index_ = ind1 * 2;
  index_2_ = ind2 * 2;
  bin_values.push_back(0);
  bin_values.push_back(0.5);
  bin_values.push_back(1.0);
  bin_values.push_back(1.5);
  bin_values.push_back(2.0);
  bin_values.push_back(2.5);
  bin_values.push_back(3.0);
  bin_values.push_back(3.5);
  bin_values.push_back(4.0);
  bin_values.push_back(4.5);
  bin_values.push_back(5.0);
  bin_values.push_back(5.5);
  bin_values.push_back(6.0);
  bin_values.push_back(6.5);
  bin_values.push_back(7.0);
  bin_values.push_back(7.5);
  bin_values.push_back(8.0);
  bin_values.push_back(8.5);
  bin_values.push_back(9.0);
  map_[0] = 0;
  map_[0.5] = 0;
  map_[1.0] = 0;
  map_[1.5] = 0;
  map_[2.0] = 0;
  map_[2.5] = 0;
  map_[3.0] = 0;
  map_[3.5] = 0;
  map_[4.0] = 0;
  map_[4.5] = 0;
  map_[5.0] = 0;
  map_[5.5] = 0;
  map_[6.0] = 0;
  map_[6.5] = 0;
  map_[7.0] = 0;
  map_[7.5] = 0;
  map_[8.0] = 0;
  map_[8.5] = 0;
  map_[9.0] = 0;
}


double
DistanceProfileRes::scale( double old_value) {
  double old_min = -1;
  double old_max = 1;
  double new_min = 180.0 * -1;
  double new_max = 180.0;

  return ( (old_value - old_min) / (old_max - old_min) ) * (new_max - new_min) + new_min;
}

void
DistanceProfileRes::add_new_calculation(const Individual& ind1,const Individual& ind2) {
  double sum = 0.0;
  int cnt = 0;
  for (int i = index_; i < index_ + 2; i++) {
    sum += std::pow(scale(ind1.vars[i]) - scale(ind2.vars[index_2_ + cnt]), 2);
    cnt++;
  }
  double euc = std::sqrt(sum);

  bool inserted = false;
  if (euc > 9.0) {
    map_[bin_values[bin_values.size() - 1]]++;
    inserted = true;
  }

  for (int i = 1; i < map_.size() && !inserted; i++) {
    if (euc < bin_values[i]) {
      map_[bin_values[i - 1]]++;
      inserted = true;
    }
  }
}

double
DistanceProfileRes::get_peak() {
  int best_bin;
  int best_count = 0;
  std::cout << "peak for " << index_  << ", " << index_2_ << std::endl;
  std::map<double, int>::iterator last = map_.end();
  --last;

  for (std::map<double, int>::iterator i = map_.begin(); i != last; i++) {
    //      std::cout << "bin " << i->first << " count " << i->second << std::endl;
    if (i->second > best_count) {
      best_count = i->second;
      best_bin = i->first;
    }
  }

  std::cout << "peak for " << index_  << ", " << index_2_ << " : is at " << best_bin << " ( " << best_count << " ) " << std::endl;

  return best_bin;
}


void
DistancesProfileCalculator::apply(const std::vector<Individual>& popul) {
  int nres = 1;
  for (int iter_res = 1; iter_res < 6; iter_res++) {
    if(iter_res != nres) {
      DistanceProfileRes distance_res(nres, iter_res);
      for (int i = 0; i < popul.size(); i++) {
	for (int j = i + 1; j < popul.size(); j++) {
	  if (i != j) {
	    distance_res.add_new_calculation(popul[i], popul[j]);
	  }
	}
      }
      double peak_res = distance_res.get_peak();
    }
  }
}
