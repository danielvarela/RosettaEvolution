
#ifndef POPULEVALUATIONTEST_H
#define POPULEVALUATIONTEST_H


#include "../Controller/DE_Operator.hpp"
#include "../MpiFiles/MasterRosettaCalculator.hpp"
#include "../Controller/StageBuilder.hpp"

class PopulEvaluationTest {
public:
  boost::shared_ptr<DE_Operator> app_operator;
  std::vector<Individual> current_population;
  boost::shared_ptr<MoverDE> de;
  boost::shared_ptr<MasterRosettaCalculator> mpi_calculator;
  boost::shared_ptr<PoseScoreFunction> simple_ffxn;

  explicit PopulEvaluationTest(boost::shared_ptr<DE_Operator> app_operator_in);

  std::vector<Individual> start_popul();

  void init_operator_for_test();

  bool run();

};




#endif /* POPULEVALUATIONTEST_H */
