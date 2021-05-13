

#ifndef CRYOPROTOCOLTEST_H
#define CRYOPROTOCOLTEST_H

#include "../Controller/DE_Operator.hpp"

class CryoProtocolTest {
public:
  boost::shared_ptr<DE_Operator> app_operator;
  std::vector<Individual> current_population;
  boost::shared_ptr<MoverDE> de;
  FitFunctionPtr scfxn;

  explicit CryoProtocolTest(boost::shared_ptr<DE_Operator> app_operator_in);

  std::vector<Individual> start_popul();

  void init_operator_for_test();

  bool run();

};





#endif /* CRYOPROTOCOLTEST_H */
