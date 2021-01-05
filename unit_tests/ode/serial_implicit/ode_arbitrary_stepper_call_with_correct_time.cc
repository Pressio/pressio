
#include <gtest/gtest.h>
#include "pressio_ode_implicit.hpp"

namespace{
struct MyApp
{
  mutable int count_ = 0;
  using scalar_type = double;
  using state_type    = Eigen::VectorXd;
  using discrete_time_residual_type = Eigen::VectorXd;
  using discrete_time_jacobian_type = Eigen::SparseMatrix<double>;

public:
  discrete_time_residual_type createDiscreteTimeResidual() const{
    discrete_time_residual_type R(3);
    return R;
  }
  discrete_time_jacobian_type createDiscreteTimeJacobian() const{
    discrete_time_jacobian_type J(3,3);
    return J;
  }

  template <typename step_t, typename ... Args>
  void discreteTimeResidual(const step_t & step,
                            const scalar_type & time,
                            const scalar_type & dt,
                            discrete_time_residual_type & R,
                            Args && ... args) const
  {
    ++count_;
    std::cout << count_ << " " << time << " " << dt << std::endl;

    EXPECT_DOUBLE_EQ(dt, 2.5);
    // count_==1: we are doing first step: t_0 -> t_1
    // so evaluation time should be at following step
    if(count_==1){
      EXPECT_EQ(step, 1);
      EXPECT_DOUBLE_EQ(time, dt);
    }

    // count_==2: we are doing step: t_1 -> t_2
    // so evaluation time should be at following step
    if(count_==2){
      EXPECT_EQ(step, 2);
      EXPECT_DOUBLE_EQ(time, 2.*dt);
    }
  }

  template <typename step_t, typename ... Args>
  void discreteTimeJacobian(const step_t & step,
                            const scalar_type & time,
                            const scalar_type & dt,
                            discrete_time_jacobian_type & J,
                            Args && ... states) const
  {
    EXPECT_DOUBLE_EQ(dt, 2.5);
    if(count_==1){
      EXPECT_EQ(step, 1);
      EXPECT_DOUBLE_EQ(time, dt);
    }
    if(count_==2){
      EXPECT_EQ(step, 2);
      EXPECT_DOUBLE_EQ(time, 2.*dt);
    }
  }
};

template<typename R_t, typename J_t>
struct MyFakeSolver
{
  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    R_t R(3);
    J_t J(3,3);
    sys.residual(state, R);
    sys.jacobian(state, J);
  }
};
}

TEST(ode_implicit_arbitrary, callWithCorrectTime1)
{
  // fake test to make sure the integrator calls
  // the model with the correct time

  using namespace pressio;
  using app_t = MyApp;
  using nstate_t	= typename app_t::state_type;
  using nresid_t	= typename app_t::discrete_time_residual_type;
  using njacobian_t	= typename app_t::discrete_time_jacobian_type;
  using state_t		= ::pressio::containers::Vector<nstate_t>;
  using res_t		= ::pressio::containers::Vector<nresid_t>;
  using jac_t		= ::pressio::containers::SparseMatrix<njacobian_t>;

  app_t appObj;
  state_t y(3); y.data()->setConstant(1.);

  using my_custom_order = ::pressio::ode::types::StepperOrder<1>;
  using my_num_states	= ::pressio::ode::types::StepperTotalNumberOfStates<2>;
  using stepper_t = ::pressio::ode::ImplicitStepper<
    ::pressio::ode::implicitmethods::Arbitrary,
    state_t, res_t, jac_t, app_t, my_custom_order, my_num_states>;
  stepper_t stepperObj(y, appObj);

  MyFakeSolver<res_t, jac_t> solver;
  double dt = 2.5;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 2, solver);
}

TEST(ode_implicit_arbitrary, callWithCorrectTime2)
{
  using namespace pressio;
  using app_t = MyApp;
  using nstate_t	= typename app_t::state_type;
  using nresid_t	= typename app_t::discrete_time_residual_type;
  using njacobian_t	= typename app_t::discrete_time_jacobian_type;
  using state_t		= ::pressio::containers::Vector<nstate_t>;
  using res_t		= ::pressio::containers::Vector<nresid_t>;
  using jac_t		= ::pressio::containers::SparseMatrix<njacobian_t>;

  app_t appObj;
  state_t y(3); y.data()->setConstant(1.);

  using my_custom_order = ::pressio::ode::types::StepperOrder<1>;
  using my_num_states	= ::pressio::ode::types::StepperTotalNumberOfStates<2>;
  using stepper_t = ::pressio::ode::ImplicitStepper<
    ::pressio::ode::implicitmethods::Arbitrary,
    state_t, res_t, jac_t, app_t, my_custom_order, my_num_states>;
  stepper_t stepperObj(y, appObj);

  MyFakeSolver<res_t, jac_t> solver;
  double dt = 2.5;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 2, solver);
}

TEST(ode_implicit_arbitrary, callWithCorrectTime3)
{
  using namespace pressio;
  using app_t = MyApp;
  using nstate_t	= typename app_t::state_type;
  using nresid_t	= typename app_t::discrete_time_residual_type;
  using njacobian_t	= typename app_t::discrete_time_jacobian_type;
  using state_t		= ::pressio::containers::Vector<nstate_t>;
  using res_t		= ::pressio::containers::Vector<nresid_t>;
  using jac_t		= ::pressio::containers::SparseMatrix<njacobian_t>;

  app_t appObj;
  state_t y(3); y.data()->setConstant(1.);

  using my_custom_order = ::pressio::ode::types::StepperOrder<1>;
  using my_num_states	= ::pressio::ode::types::StepperTotalNumberOfStates<3>;
  using stepper_t = ::pressio::ode::ImplicitStepper<
    ::pressio::ode::implicitmethods::Arbitrary,
    state_t, res_t, jac_t, app_t, my_custom_order, my_num_states>;
  stepper_t stepperObj(y, appObj);

  MyFakeSolver<res_t, jac_t> solver;
  double dt = 2.5;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 2, solver);
}
