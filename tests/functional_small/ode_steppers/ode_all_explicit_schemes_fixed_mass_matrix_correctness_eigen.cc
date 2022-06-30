
#include <gtest/gtest.h>
#include "pressio/ode_steppers_explicit.hpp"
#include "pressio/ode_advancers.hpp"

struct MyApp2WithMM
{
  using independent_variable_type = double;
  using state_type           = Eigen::VectorXd;
  using right_hand_side_type = state_type;
  using mass_matrix_type     = Eigen::MatrixXd;

  mutable int count1 = 0;
  std::map<int, Eigen::VectorXd> & rhs_;
  const Eigen::MatrixXd & uniqueMM_;

  MyApp2WithMM(std::map<int, Eigen::VectorXd> & rhs,
	       const Eigen::MatrixXd & MM)
    : rhs_(rhs), uniqueMM_(MM){}

  state_type createState() const{
    state_type ret(3); ret.setZero();
    return ret;
  }

  right_hand_side_type createRightHandSide() const{
    right_hand_side_type ret(3); ret.setZero();
    return ret;
  };

  mass_matrix_type createMassMatrix() const{
    mass_matrix_type ret(3,3); ret.setZero();
    return ret;
  };

  void rightHandSide(const state_type & y,
		     independent_variable_type evaltime,
		     right_hand_side_type & rhs) const
  {
    rhs = Eigen::VectorXd::Random(rhs.rows());
    for (int i=0; i<rhs.size(); ++i){
      rhs(i) = evaltime;
    }

    rhs_[++count1] = rhs;
  };

  void massMatrix(mass_matrix_type & M) const
  {
    M = uniqueMM_;
  };
};

struct MyApp2NoMM
{
  using independent_variable_type = double;
  using state_type           = Eigen::VectorXd;
  using right_hand_side_type = state_type;

  mutable int count1 = 0;
  const std::map<int, Eigen::VectorXd> & rhs_;
  const Eigen::MatrixXd & uniqueMM_;

  MyApp2NoMM(const std::map<int, Eigen::VectorXd> & rhs,
	       const Eigen::MatrixXd & MM)
    : rhs_(rhs), uniqueMM_(MM){}

  state_type createState() const{
    state_type ret(3); ret.setZero();
    return ret;
  }

  right_hand_side_type createRightHandSide() const{
    right_hand_side_type ret(3); ret.setZero();
    return ret;
  };

  void rightHandSide(const state_type & y,
		     independent_variable_type evaltime,
		     right_hand_side_type & rhs) const
  {
    rhs = rhs_.at(++count1);
    // we want to make sure we have the correct time
    for (int i=0; i<rhs.size(); ++i){
      EXPECT_DOUBLE_EQ(rhs(i), evaltime);
    }

    rhs = uniqueMM_.inverse()*rhs;
  };
};

struct LinearSolver2
{
  void solve(const Eigen::MatrixXd & A,
	     Eigen::VectorXd & x,
	     const Eigen::VectorXd & b)
  {
    x = A.colPivHouseholderQr().solve(b);
  }
};

#define ODE_MASS_MATRIX_CHECK_TEST(NAME)			\
  /* for given method, M dy/dt = f with mass-matrix API */		\
  /* should be same as doing dy/dt = M^-1 f */				\
  std::cout << "\n";							\
  using namespace pressio;						\
  srand(342556331);							\
  const auto nsteps = ::pressio::ode::StepCount(4);			\
  const double dt = 2.;							\
  std::map<int, Eigen::VectorXd> rhs;					\
  /**/									\
  /* fix MM once */							\
  Eigen::MatrixXd M = Eigen::MatrixXd::Random(3, 3);			\
  /**/									\
  /* first, solve problem using mass matrix API */			\
  Eigen::VectorXd y0(3);						\
  y0.setZero();{							\
    MyApp2WithMM appObj(rhs, M);					\
    auto stepperObj = ode::create_##NAME##_stepper(appObj);		\
    LinearSolver2 solver;						\
    y0(0) = 1.; y0(1) = 2.; y0(2) = 3.;					\
    ode::advance_n_steps(stepperObj, y0, 0.0, dt, nsteps, solver);	\
  }									\
  std::cout << "y0 : \n" << y0 << " \n";				\
  /**/									\
  /* solve problem computing modified rhs using mass matrix inverse */	\
  Eigen::VectorXd y1(3);						\
  y1.setZero();{							\
    MyApp2NoMM appObj(rhs, M);						\
    auto stepperObj = ode::create_##NAME##_stepper(appObj);		\
    y1(0) = 1.; y1(1) = 2.; y1(2) = 3.;					\
    ode::advance_n_steps(stepperObj, y1, 0.0, dt, nsteps);		\
  }									\
  std::cout << "y1 : \n" << y1 << " \n";				\
  /**/									\
  EXPECT_NEAR( y0(0), y1(0), 1e-12);					\
  EXPECT_NEAR( y0(1), y1(1), 1e-12);					\
  EXPECT_NEAR( y0(2), y1(2), 1e-12);					\


TEST(ode_explicit_steppers, forward_euler_with_mass_matrix_use_inverse){
  ODE_MASS_MATRIX_CHECK_TEST(forward_euler)
}

TEST(ode_explicit_steppers, ab2_with_mass_matrix_use_inverse){
  ODE_MASS_MATRIX_CHECK_TEST(ab2)
}

TEST(ode_explicit_steppers, rk4_with_mass_matrix_use_inverse){
  ODE_MASS_MATRIX_CHECK_TEST(rk4)
}

TEST(ode_explicit_steppers, spprk3_with_mass_matrix_use_inverse){
  ODE_MASS_MATRIX_CHECK_TEST(ssprk3)
}
