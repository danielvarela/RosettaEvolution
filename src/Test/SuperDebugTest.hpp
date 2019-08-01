
#ifndef SUPERDEBUGTEST_H
#define SUPERDEBUGTEST_H

#include "../Controller/DE_Operator.hpp"
#include "../MpiFiles/MasterRosettaCalculator.hpp"
#include "../Controller/StageBuilder.hpp"

class SuperDebugTest {
public:
  boost::shared_ptr<DE_Operator> app_operator;
  std::vector<Individual> current_population;
  boost::shared_ptr<MoverDE> de;
  boost::shared_ptr<MasterRosettaCalculator> mpi_calculator;
  boost::shared_ptr<PoseScoreFunction> simple_ffxn;
  core::pose::PoseOP native_pose;

  explicit SuperDebugTest(boost::shared_ptr<DE_Operator> app_operator_in);

  std::vector<Individual> start_popul();

  void init_operator_for_test();

  bool run();

  void make_population_straight();
  void insert_native_pose_at_popul();

};




#endif /* SUPERDEBUGTEST_H */
