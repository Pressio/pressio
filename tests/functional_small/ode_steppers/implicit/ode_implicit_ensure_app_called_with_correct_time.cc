
#include <gtest/gtest.h>
#include "pressio/ode_steppers_implicit.hpp"
#include "pressio/ode_advancers.hpp"

namespace
{
struct MyApp1
{
  using independent_variable_type   = double;
  using state_type    = Eigen::VectorXd;
  using right_hand_side_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

public:
  state_type createState() const{
    return state_type(3);
  }

  right_hand_side_type createRightHandSide() const{
    right_hand_side_type f(3);
    return f;
  }

  jacobian_type createJacobian() const{
    jacobian_type J(3,3);
    return J;
  }

  void rightHandSide(const state_type & yn,
		const independent_variable_type& evaltime,
		right_hand_side_type & f) const
  {
    std::cout << "f: t=" << evaltime << "\n";
    f[0] = evaltime+1.;
    f[1] = evaltime+2.;
    f[2] = evaltime+3.;
  }

  void jacobian(const state_type&, const independent_variable_type&, jacobian_type&) const{
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
  app_t appObj;
  state_t y(3); y.setConstant(1.5);

  auto stepperObj = ode::create_bdf1_stepper(appObj);
  MyFakeSolver1 solver;
  double dt = 1.1;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, ::pressio::ode::StepCount(1), solver);
}


namespace
{
struct MyApp2
{
  using independent_variable_type   = double;
  using state_type    = Eigen::VectorXd;
  using right_hand_side_type = state_type;
  using jacobian_type = Eigen::SparseMatrix<double>;

public:
  state_type createState() const{
    return state_type(3);
  }

  right_hand_side_type createRightHandSide() const{
    right_hand_side_type f(3);
    return f;
  }

  jacobian_type createJacobian() const{
    jacobian_type J(3,3);
    return J;
  }

  void rightHandSide(const state_type & yn,
    const independent_variable_type& evaltime,
    right_hand_side_type & f) const
  {
    std::cout << "f: t=" << evaltime << "\n";
    f[0] = evaltime+1.;
    f[1] = evaltime+2.;
    f[2] = evaltime+3.;
  }

  void jacobian(const state_type&, const independent_variable_type&, jacobian_type&) const{
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
  app_t appObj;
  state_t y(3); y.setConstant(1.5);

  auto stepperObj = ode::create_bdf2_stepper(appObj);
  MyFakeSolver2 solver;
  double dt = 1.1;
  ode::advance_n_steps(stepperObj, y, 0.0, dt, ::pressio::ode::StepCount(2), solver);
}


// namespace
// {
// struct MyApp3
// {
//   mutable int count_ = 0;
//   using independent_variable_type = double;
//   using state_type    = Eigen::VectorXd;
//   using discrete_time_residual_type = Eigen::VectorXd;
//   using discrete_time_jacobian_type = Eigen::SparseMatrix<double>;

// public:
//   state_type createState() const{
//     return state_type(3);
//   }

//   discrete_time_residual_type createDiscreteTimeResidual() const{
//     discrete_time_residual_type R(3);
//     return R;
//   }
//   discrete_time_jacobian_type createDiscreteTimeJacobian() const{
//     discrete_time_jacobian_type J(3,3);
//     return J;
//   }

//   template <typename step_t>
//   void discreteTimeResidual(const step_t & step,
//                             const independent_variable_type & time,
//                             const independent_variable_type & dt,
//                             discrete_time_residual_type & R,
//                             const state_type &,
//                             const state_type &) const
//   {
//     ++count_;
//     std::cout << count_ << " " << time << " " << dt << std::endl;

//     EXPECT_DOUBLE_EQ(dt, 2.5);
//     // count_==1: we are doing first step: t_0 -> t_1
//     // so evaluation time should be at following step
//     if(count_==1){
//       EXPECT_EQ(step, 1);
//       EXPECT_DOUBLE_EQ(time, dt);
//     }

//     // count_==2: we are doing step: t_1 -> t_2
//     // so evaluation time should be at following step
//     if(count_==2){
//       EXPECT_EQ(step, 2);
//       EXPECT_DOUBLE_EQ(time, 2.*dt);
//     }
//   }

//   template <typename step_t>
//   void discreteTimeJacobian(const step_t & step,
//                             const independent_variable_type & time,
//                             const independent_variable_type & dt,
//                             discrete_time_jacobian_type & J,
//                             const state_type &,
//                             const state_type &) const
//   {
//     EXPECT_DOUBLE_EQ(dt, 2.5);
//     if(count_==1){
//       EXPECT_EQ(step, 1);
//       EXPECT_DOUBLE_EQ(time, dt);
//     }
//     if(count_==2){
//       EXPECT_EQ(step, 2);
//       EXPECT_DOUBLE_EQ(time, 2.*dt);
//     }
//   }
// };

// template<typename R_t, typename J_t>
// struct MyFakeSolver3
// {
//   template<typename system_t, typename state_t>
//   void solve(const system_t & sys, state_t & state)
//   {
//     R_t R(3);
//     J_t J(3,3);
//     sys.residual(state, R);
//     sys.jacobian(state, J);
//   }
// };
// }

// TEST(ode, mplicit_arbitrary_callWithCorrectTime1)
// {
//   // fake test to make sure the integrator calls
//   // the model with the correct time

//   using namespace pressio;
//   using app_t = MyApp3;
//   using state_t  = typename app_t::state_type;
//   app_t appObj;
//   state_t y(3); y.setConstant(1.);

//   auto stepperObj = ode::create_arbitrary_stepper<2>(appObj);

//   using res_t  = typename app_t::discrete_time_residual_type;
//   using jac_t = typename app_t::discrete_time_jacobian_type;
//   MyFakeSolver3<res_t, jac_t> solver;
//   double dt = 2.5;
//   ode::advance_n_steps(stepperObj, y, 0.0, dt, ::pressio::ode::StepCount(2), solver);
// }

// TEST(ode, implicit_arbitrary_callWithCorrectTime2)
// {
//   using namespace pressio;
//   using app_t = MyApp3;
//   using state_t  = typename app_t::state_type;
//   using res_t  = typename app_t::discrete_time_residual_type;
//   using jac_t = typename app_t::discrete_time_jacobian_type;
//   app_t appObj;
//   state_t y(3); y.setConstant(1.);

//   auto stepperObj = ode::create_arbitrary_stepper<2>(appObj);
//   MyFakeSolver3<res_t, jac_t> solver;
//   double dt = 2.5;
//   ode::advance_n_steps(stepperObj, y, 0.0, dt, ::pressio::ode::StepCount(2), solver);
// }

// TEST(ode, implicit_arbitrary_callWithCorrectTime3)
// {
//   using namespace pressio;
//   using app_t = MyApp3;
//   using state_t  = typename app_t::state_type;
//   using res_t  = typename app_t::discrete_time_residual_type;
//   using jac_t = typename app_t::discrete_time_jacobian_type;
//   app_t appObj;
//   state_t y(3); y.setConstant(1.);

//   auto stepperObj = ode::create_arbitrary_stepper<2>(appObj);
//   MyFakeSolver3<res_t, jac_t> solver;
//   double dt = 2.5;
//   ode::advance_n_steps(stepperObj, y, 0.0, dt, ::pressio::ode::StepCount(2), solver);
// }
