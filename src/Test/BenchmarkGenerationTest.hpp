
#ifndef BENCHMARKGENERATIONTEST_H
#define BENCHMARKGENERATIONTEST_H

#include "../Controller/DE_Operator.hpp"
#include "../Controller/StageBuilder.hpp"

class BenchmarkGenerationTest {
public:
  boost::shared_ptr<DE_Operator> app_operator;
  std::vector<Individual> current_population;
  boost::shared_ptr<MoverDE> de;

  explicit BenchmarkGenerationTest(boost::shared_ptr<DE_Operator> app_operator_in);

  std::vector<Individual> start_popul();

  void init_operator_for_test();

  bool run();

};




#endif /* BENCHMARKGENERATIONTEST_H */
