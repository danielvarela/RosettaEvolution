#ifndef DIFFERENTIALEVOLUTION_FITFUNCTION_HPP
#define DIFFERENTIALEVOLUTION_FITFUNCTION_HPP

#include "../Algorithm/DE_types.hpp"
#include <cmath>
#include <string>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <boost/shared_ptr.hpp>

template<class K, class V>
class MyFuncStatsMap : public std::map<K, V>
{
public:
  std::map<std::string, double> stats;

  MyFuncStatsMap() {}

  double& operator[](std::string index) {
    return stats[index];
  }

  void clear() {
    stats.clear();
  }

  std::map<std::string, double>::iterator find( const std::string& key ) {
    return stats.find(key);
  }


  void increment(std::string stat_name) {
    if (stats.find(stat_name) != stats.end()) {
	stats[stat_name]++;
      } else {
	stats[stat_name] = 1;
      }
  }

  void record(std::string stat_name, double value) {
  if (stats.find(stat_name) != stats.end()) {
	stats[stat_name] += value;
      } else {
	stats[stat_name] = value;
      }
  }

  double obtain(std::string stat_name) {
    double value;
    if (stats.find(stat_name) != stats.end()) {
      value = stats[stat_name];
    } else {
      value = 0;
    }
    return value;
  }
};

typedef MyFuncStatsMap<std::string, double> FuncStats;

class FitFunction
{
public:
  FuncStats stats;
  std::string name_;

  FitFunction(){
    old_min = -1;
    old_max = 1;
  }

  virtual double score(Individual& ind) = 0;
  virtual int D() = 0;
  void set_name(std::string input) {
    name_ = input;
  }
  virtual std::string name() = 0;

  virtual double lim() = 0;

  virtual void fill_pose(core::pose::PoseOP& p, const Individual& ind, std::string ss_) {
  }

  double scale( double old_value) ;

  virtual Individual convert(const Individual& ind);

  void reset_stats() {
    stats.clear();
  }

  virtual std::string print_stats() {
    return "no_stats";
  }

  virtual int start_res() {
    return 1;
  }

protected:
  double old_min;
  double old_max;
};

typedef boost::shared_ptr<FitFunction> FitFunctionPtr;

class DeJongFunction : public FitFunction
{
public:
  DeJongFunction() : FitFunction() {
    name_ = std::string("first De Jong");
  }

  double score(Individual& ind) {
    Individual ind_converted = convert(ind);

    double result = 0.0;
    for (int i =0 ; i < D(); ++i) {
      result += std::pow(ind_converted.vars[i], 2);
    }

    ind.score = result;
    return result;
  }

  int D() {
    return 2;
  }

  double lim() {
    return 5.12;
  }


  std::string name() {
    return name_;
  }

};

class HimmelblauFunction : public FitFunction
{
public:
  double score(Individual& ind) {
    Individual ind_converted = convert(ind);

    double result = 0.0;
    result = std::pow( std::pow(ind_converted.vars[0], 2) + ind_converted.vars[1] - 11 , 2) +
      std::pow( ind_converted.vars[0] + std::pow( ind_converted.vars[1], 2) - 7 , 2) ;

    ind.score = result;
    return result;
  }

  int D() {
    return 2;
  }

  double lim() {
    return 6;
  }


  std::string name() {
    return std::string("first De Jong");
  }
};


class ToyFunction : public FitFunction
{
public:
  double score(Individual& ind) {
    Individual ind_converted = convert(ind);

    double result = 0.0;
    result = -0.5 + std::pow(ind_converted.vars[0], 2);
    ind.score = result;
    return result;
  }

  int D() {
    return 1;
  }

  double lim() {
    return 1;
  }


  std::string name() {
    return std::string("first De Jong");
  }
};



#endif // DIFFERENTIALEVOLUTION_FITFUNCTION_HPP
