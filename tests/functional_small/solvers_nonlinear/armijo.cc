#include <gtest/gtest.h>
#include "pressio/solvers_nonlinear.hpp"

using _this_state_type = Eigen::VectorXd;
constexpr int _this_N = 3;

/*
 "simple test to make sure armijo is working correct but without any meaning
*/
struct FakeSolver
{
  mutable int count_ = 0;
  _this_state_type correction_;
  _this_state_type gradient_;
  double residNormCurrCorrStep_ = {};
  int testCaseId_ = 0;

  FakeSolver(int testCaseId)
    : correction_(_this_N),
      gradient_(_this_N) ,
      testCaseId_(testCaseId)
  {
    correction_(0) = 1.;
    correction_(1) = 2.;
    correction_(2) = 3.;
    gradient_(0) = 1.;
    gradient_(1) = 1.;
    gradient_(2) = 1.;
  }

  int count() const { return count_; }

  const _this_state_type & correctionCRef() const{ return correction_; }
  const _this_state_type & gradientCRef() const{ return gradient_; }

  double residualNormCurrentCorrectionStep() const{ return 2.5; }

  template<class SystemType>
  void residualNorm(const SystemType & systemObj,
		    const _this_state_type & state,
		    double & residualNorm) const
  {
    if (testCaseId_ == 1){
      if (count_ == 0){ residualNorm = 3.0; }
      else if (count_ == 1){ residualNorm = 5.0; }
      else if (count_ == 2){ residualNorm = 8.0; }
      else{ residualNorm = 1.0; }
    }
    else if (testCaseId_ == 2){
      residualNorm = 11.0;
    }
    count_++;
  }
};

struct DummySystem{};

TEST(solvers_nonlinear, armijo_test1)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  _this_state_type state(_this_N);
  state.setConstant(5.);
  DummySystem sys;

  pressio::nonlinearsolvers::impl::ArmijoUpdater<_this_state_type> updater(state);
  FakeSolver solver(1);
  updater(sys, state, solver);
  std::cout << state << std::endl;

  EXPECT_TRUE(solver.count() == 4);
  EXPECT_TRUE(state(0) == 5. + 0.125 * 1.);
  EXPECT_TRUE(state(1) == 5. + 0.125 * 2.);
  EXPECT_TRUE(state(2) == 5. + 0.125 * 3.);

  pressio::log::finalize();
}

TEST(solvers_nonlinear, armijo_test2)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  _this_state_type state(_this_N);
  state.setConstant(5.);
  DummySystem sys;

  pressio::nonlinearsolvers::impl::ArmijoUpdater<_this_state_type> updater(state);
  FakeSolver solver(2);
  try{
    updater(sys, state, solver);
  }
  catch (::pressio::eh::LineSearchStepTooSmall const &e) {
    PRESSIOLOG_WARN(e.what());
  }
  std::cout << state << std::endl;

  EXPECT_TRUE(solver.count() == 10);

  // solution should not change because this test intentionally makes the update fail
  EXPECT_TRUE(state(0) == 5.);
  EXPECT_TRUE(state(1) == 5.);
  EXPECT_TRUE(state(2) == 5.);

  pressio::log::finalize();
}
