
#include <gtest/gtest.h>
#include "pressio_ode_implicit.hpp"

namespace
{
struct MyApp1
{
  using scalar_type   = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

public:
  velocity_type createVelocity() const{
    velocity_type f(3);
    return f;
  }

  jacobian_type createJacobian() const{
    jacobian_type J(3,3);
    return J;
  }

  void velocity(const state_type & yn,
		const scalar_type& time,
		velocity_type & f) const
  {
    std::cout << "f: t=" << time << "\n";
    f[0] = time+1.;
    f[1] = time+2.;
    f[2] = time+3.;
  }

  void jacobian(const state_type&, const scalar_type&, jacobian_type&) const{
    // not used by the test
  }
};

struct MyFakeSolver1
{
  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    state_t R(3);
    sys.residual(state, R);

    // here we should have:
    // R = y_n - y_n-1 - dt*f()
    // computedd at t_1 = dt = 1.1
    // where y_n = y_n-1 = [1.5 1.5 1.5]
    // f = [2.1 3.1 4.1]
    //
    // so R = -1.1*[2.1 3.1 4.1] = -[2.31 3.41 4.51]
    std::cout << R << std::endl;
    EXPECT_DOUBLE_EQ(R(0), -2.31);
    EXPECT_DOUBLE_EQ(R(1), -3.41);
    EXPECT_DOUBLE_EQ(R(2), -4.51);
  }
};
}

TEST(ode, implicit_euler_appVelocityCalledWithCorrectTime)
{
  // fake test to make sure the app velocity
  // for bdf1 is called with the time at the next time step
  // so here we set the velocity to be time+something
  // so that when the "fake" solver is called we can verify
  // that the residual has the correct numbers

  using namespace pressio;
  using app_t = MyApp1;
  using state_t = typename app_t::state_type;
  using res_t = typename app_t::velocity_type;
  using jac_t = typename app_t::jacobian_type;
  app_t appObj;
  state_t y(3); y.setConstant(1.5);

  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::Euler,
    state_t, res_t, jac_t, app_t>;
  stepper_t stepperObj(y, appObj);

  MyFakeSolver1 solver;
  double dt = 1.1;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 1, solver);
}


namespace
{
struct MyApp2
{
  using scalar_type   = double;
  using state_type    = Eigen::VectorXd;
  using velocity_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

public:
  velocity_type createVelocity() const{
    velocity_type f(3);
    return f;
  }

  jacobian_type createJacobian() const{
    jacobian_type J(3,3);
    return J;
  }

  void velocity(const state_type & yn,
    const scalar_type& time,
    velocity_type & f) const
  {
    std::cout << "f: t=" << time << "\n";
    f[0] = time+1.;
    f[1] = time+2.;
    f[2] = time+3.;
  }

  void jacobian(const state_type&, const scalar_type&, jacobian_type&) const{
    // not used by the test
  }
};

struct MyFakeSolver2
{
  int count_ = {};

  template<typename system_t, typename state_t>
  void solve(const system_t & sys, state_t & state)
  {
    ++count_;

    state_t R(3);
    sys.residual(state, R);

    if (count_==1){
      // this is called from auxiliary stepper which is bdf1

      // here we should have:
      // R = y_n - y_n-1 - dt*f()
      // computedd at t_1 = dt = 1.1
      // where y_n = y_n-1 = [1.5 1.5 1.5]
      // f = [2.1 3.1 4.1]
      EXPECT_DOUBLE_EQ(R(0), 1.5-1.5-1.1*2.1);
      EXPECT_DOUBLE_EQ(R(1), 1.5-1.5-1.1*3.1);
      EXPECT_DOUBLE_EQ(R(2), 1.5-1.5-1.1*4.1);
    }
    if (count_==2)
    {
      // this is called from bdf2
      // here we should have:
      // R = y_n+1 - (4/3)y_n + (1/3)*y_n-1 - (2/3)*dt*f()
      // computedd at t_2 = dt = 2.2
      // where y_n+1 = y_n = y_n-1 = [1.5 1.5 1.5]
      // f = [3.2 4.2 5.2]
      //
      std::cout << state << std::endl;
      std::cout << R << std::endl;
      EXPECT_DOUBLE_EQ(R(0), 1.5-(4./3.)*1.5+(1./3.)*1.5 - (2./3.)*1.1*3.2);
      EXPECT_DOUBLE_EQ(R(1), 1.5-(4./3.)*1.5+(1./3.)*1.5 - (2./3.)*1.1*4.2);
      EXPECT_DOUBLE_EQ(R(2), 1.5-(4./3.)*1.5+(1./3.)*1.5 - (2./3.)*1.1*5.2);
    }
  }
};
}

TEST(ode, implicit_bdf2_appVelocityCalledWithCorrectTime)
{
  // fake test to make sure the app velocity
  // for bdf1 is called with the time at the next time step
  // so here we set the velocity to be time+something
  // so that when the "fake" solver is called we can verify
  // that the residual has the correct numbers

  using namespace pressio;
  using app_t = MyApp2;
  using state_t = typename app_t::state_type;
  using res_t = typename app_t::velocity_type;
  using jac_t = typename app_t::jacobian_type;
  app_t appObj;
  state_t y(3); y.setConstant(1.5);

  using aux_stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::Euler, state_t, res_t, jac_t, app_t>;
  aux_stepper_t stepperAux(y, appObj);

  using stepper_t = ode::ImplicitStepper<
    ode::implicitmethods::BDF2,
    state_t, res_t, jac_t, app_t, aux_stepper_t>;
  stepper_t stepperObj(y, appObj, stepperAux);

  MyFakeSolver2 solver;
  double dt = 1.1;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 2, solver);
}


namespace
{
struct MyApp3
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
struct MyFakeSolver3
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

TEST(ode, mplicit_arbitrary_callWithCorrectTime1)
{
  // fake test to make sure the integrator calls
  // the model with the correct time

  using namespace pressio;
  using app_t = MyApp3;
  using state_t  = typename app_t::state_type;
  using res_t  = typename app_t::discrete_time_residual_type;
  using jac_t = typename app_t::discrete_time_jacobian_type;
  app_t appObj;
  state_t y(3); y.setConstant(1.);

  using my_custom_order = ::pressio::ode::StepperOrder<1>;
  using my_num_states = ::pressio::ode::StepperTotalNumberOfStates<2>;
  using stepper_t = ::pressio::ode::ImplicitStepper<
    ::pressio::ode::implicitmethods::Arbitrary,
    state_t, res_t, jac_t, app_t, my_custom_order, my_num_states>;
  stepper_t stepperObj(y, appObj);

  MyFakeSolver3<res_t, jac_t> solver;
  double dt = 2.5;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 2, solver);
}

TEST(ode, implicit_arbitrary_callWithCorrectTime2)
{
  using namespace pressio;
  using app_t = MyApp3;
  using state_t  = typename app_t::state_type;
  using res_t  = typename app_t::discrete_time_residual_type;
  using jac_t = typename app_t::discrete_time_jacobian_type;
  app_t appObj;
  state_t y(3); y.setConstant(1.);

  using my_custom_order = ::pressio::ode::StepperOrder<1>;
  using my_num_states = ::pressio::ode::StepperTotalNumberOfStates<2>;
  using stepper_t = ::pressio::ode::ImplicitStepper<
    ::pressio::ode::implicitmethods::Arbitrary,
    state_t, res_t, jac_t, app_t, my_custom_order, my_num_states>;
  stepper_t stepperObj(y, appObj);

  MyFakeSolver3<res_t, jac_t> solver;
  double dt = 2.5;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 2, solver);
}

TEST(ode, implicit_arbitrary_callWithCorrectTime3)
{
  using namespace pressio;
  using app_t = MyApp3;
  using state_t  = typename app_t::state_type;
  using res_t  = typename app_t::discrete_time_residual_type;
  using jac_t = typename app_t::discrete_time_jacobian_type;
  app_t appObj;
  state_t y(3); y.setConstant(1.);

  using my_custom_order = ::pressio::ode::StepperOrder<1>;
  using my_num_states = ::pressio::ode::StepperTotalNumberOfStates<3>;
  using stepper_t = ::pressio::ode::ImplicitStepper<
    ::pressio::ode::implicitmethods::Arbitrary,
    state_t, res_t, jac_t, app_t, my_custom_order, my_num_states>;
  stepper_t stepperObj(y, appObj);

  MyFakeSolver3<res_t, jac_t> solver;
  double dt = 2.5;
  ode::advanceNSteps(stepperObj, y, 0.0, dt, 2, solver);
}