

#ifndef DE_TYPES
#define DE_TYPES

#include <stdlib.h>
#include <string>
#include <vector>

static double DE_LIMIT = 1.0;
static double SCORE_ERROR_FIXED = -1000;

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

  IndImproved() {
    score = 1000;
  }

  explicit IndImproved(int size_) {
    vars.resize(size_);
    score = 1000;
  }


  IndImproved(int size_, std::string ss_) {
    vars.resize(size_ + 1);
    score = 1000;
    ss = ss_;
  }

};

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
