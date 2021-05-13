

#ifndef DE_TYPES
#define DE_TYPES

#include <stdlib.h>
#include <string>
#include <vector>

#if(MPI_ENABLED)
#include <boost/mpi.hpp>
#endif

static double DE_LIMIT = 1.0;
static double SCORE_ERROR_FIXED = -10000;

inline double URAND() {
  return (double)(rand()/(static_cast<double>(RAND_MAX)));
}

struct Parents
{
  int x1, x2, x3;

  Parents() : x1(0), x2(0), x3(0) {
  }
};

class IndImproved {
public:
  std::vector<double> vars;
  std::vector<double> omega;
  std::string ss;
  double score;
  int gen_acc;

  IndImproved() {
    score = 1000;
    gen_acc = 0;
  }

  explicit IndImproved(int size_) {
    vars.resize(size_);
    score = 1000;
    gen_acc = 0;
  }


  IndImproved(int size_, std::string ss_) {
    vars.resize(size_ + 1);
    score = 1000;
    ss = ss_;
    gen_acc = 0;
  }

#if(MPI_ENABLED)
private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar & vars;
    ar & omega;
    ar & ss;
    ar & score;
    ar & gen_acc;
  }

#endif
};

#if(MPI_ENABLED)
class IndMPI {
public:
  IndImproved ind;
  int index;
  int nearest;

  IndMPI() {}
  IndMPI(int i, IndImproved ind_in): index(i), ind(ind_in) {}
  IndMPI(int i, IndImproved ind_in, int near): index(i), ind(ind_in), nearest(near) {}

private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar & ind;
    ar & index;
    ar & nearest;
  }
};

#endif

class IndPose {
public:
  std::vector<double> vars;
  std::vector<double> omega;
  std::string ss;
  double score;

  IndPose() {
    score = 1000;
  }

  explicit IndPose(int size_) {
    vars.resize(size_);
    score = 1000;
  }


  IndPose(int size_, std::string ss_) {
    vars.resize(size_ + 1);
    score = 1000;
    ss = ss_;
  }

};


typedef IndImproved Individual;




#endif // DE_TYPES
