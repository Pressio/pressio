#include <gtest/gtest.h>
#include "pressio/solvers_nonlinear.hpp"
using _this_state_type = Eigen::VectorXd;
constexpr int _this_N = 3;

struct FakeSolver
{
  mutable int count_ = 0;
  _this_state_type correction_;
  double residNormCurrCorrStep_ = {};

  FakeSolver() : correction_(_this_N){
    correction_(0) = 1.;
    correction_(1) = 2.;
    correction_(2) = 3.;
  }

  int count() const { return count_; }

  const _this_state_type & correctionCRef() const{
    return correction_;
  }

  double residualNormCurrentCorrectionStep() const{ return 2.5; }

  template<class SystemType>
  void residualNorm(const SystemType & systemObj,
		    const _this_state_type & state,
		    double & residualNorm) const
  {
    if (count_ == 0){
      residualNorm = 3.0;
    }
    else if (count_ == 1){
      residualNorm = 5.0;
    }
    else{
      residualNorm = 1.0;
    }
    count_++;
  }

};

struct DummySystem{};

TEST(solvers_nonlinear, backtracking)
{
  pressio::log::initialize(pressio::logto::terminal);
  pressio::log::setVerbosity({pressio::log::level::debug});

  _this_state_type state(_this_N);
  state.setConstant(5.);
  DummySystem sys;

  pressio::nonlinearsolvers::impl::BacktrackStrictlyDecreasingObjectiveUpdater<_this_state_type> updater(state);
  FakeSolver solver;
  updater(sys, state, solver);

  EXPECT_TRUE(solver.count() == 3);
  EXPECT_TRUE(state(0) == 5. + 0.25 * 1.);
  EXPECT_TRUE(state(1) == 5. + 0.25 * 2.);
  EXPECT_TRUE(state(2) == 5. + 0.25 * 3.);

  pressio::log::finalize();
}
