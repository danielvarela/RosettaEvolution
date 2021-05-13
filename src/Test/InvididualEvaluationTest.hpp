#ifndef INVIDIDUALEVALUATIONTEST_H
#define INVIDIDUALEVALUATIONTEST_H


#include "../Controller/DE_Operator.hpp"

class IndividualEvaluationTest {
public:
  boost::shared_ptr<DE_Operator> app_operator;
  std::vector<Individual> current_population;
  boost::shared_ptr<MoverDE> de;
  FitFunctionPtr scfxn;

  explicit IndividualEvaluationTest(boost::shared_ptr<DE_Operator> app_operator_in);

  std::vector<Individual> start_popul();

  void init_operator_for_test();

  bool run();

};



#endif /* INVIDIDUALEVALUATIONTEST_H */
