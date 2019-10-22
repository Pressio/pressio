
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "SOLVERS_NONLINEAR"
#include "../reference_apps_for_testing.hpp"

template<typename state_type, typename system_type, typename residual_type>
class ResidualPolicy
  : public ::pressio::ode::policy::ImplicitResidualPolicyBase<
  ResidualPolicy<state_type, system_type, residual_type>
  >
{
public:
  static constexpr auto stepper_order = 1;
  static constexpr auto num_aux_states = 1;

  void operator()(const state_type & y,
		  residual_type & R,
		  const std::array<state_type, num_aux_states> & oldYs,
		  const system_type & model,
		  double t,
		  double dt,
		  ::pressio::ode::types::step_t step) const
  {
    model.velocity(*y.data(), t, *R.data());
    R.data()->setConstant(1);
  }

  residual_type operator()(const state_type & y,
  			   const std::array<state_type, num_aux_states> & oldYs,
  			   const system_type & model,
  			   double t,
  			   double dt,
			   ::pressio::ode::types::step_t step) const{
    // here I would need to compute the time discrete residual
    residual_type R(3);
    R.data()->setConstant(1);
    return R;
  }
};//end class


template<typename state_type, typename system_type, typename jacobian_type>
class JacobianPolicy
  : public ::pressio::ode::policy::JacobianPolicyBase<
  JacobianPolicy<state_type, system_type, jacobian_type>
  >
{
public:
  static constexpr auto stepper_order = 1;
  static constexpr auto num_aux_states = 1;

  void operator()(const state_type & y,
		  jacobian_type & J,
		  const system_type & model,
		  double t,
		  double dt,
		  ::pressio::ode::types::step_t step) const{
    J.resize(3,3);
    // here I would need to compute the time discrete version
  }

  jacobian_type operator()(const state_type & y,
  			   const system_type & model,
  			   double t,
  			   double dt,
			   ::pressio::ode::types::step_t step) const{
    // here I would need to compute the time discrete version
    jacobian_type J(3, 3);
    return J;
  }

};//end class


TEST(ode_implicit, arbitraryStepperRunEuler)
{
  using namespace pressio;
  using app_t = ode::testing::refAppForImpEigen;
  using nstate_t = typename app_t::state_type;
  using nveloc_t = typename app_t::velocity_type;
  using njacobian_t = typename app_t::jacobian_type;
  app_t appObj;

  using state_t = containers::Vector<nstate_t>;
  using res_t = containers::Vector<nveloc_t>;
  using jac_t = containers::Matrix<njacobian_t>;

  using residual_policy_t = ResidualPolicy<state_t, app_t, res_t>;
  using jacobian_policy_t = JacobianPolicy<state_t, app_t, jac_t>;
  static_assert
    (ode::meta::is_legitimate_residual_policy_for_implicit_arbitrary_stepper<
     residual_policy_t, state_t, res_t, app_t, double>::value, "");
  static_assert
    (ode::meta::is_legitimate_jacobian_policy_for_implicit_arbitrary_stepper<
     jacobian_policy_t, state_t, jac_t, app_t, double>::value, "");

  using stepper_t = ode::ImplicitStepper<
    ode::ImplicitEnum::Arbitrary,
    state_t, res_t, jac_t, app_t, residual_policy_t, jacobian_policy_t>;


  state_t y(3);
  y[0] = 1.; y[1] = 2.; y[2] = 3.;

  residual_policy_t resPolObj;
  jacobian_policy_t jacPolObj;
  stepper_t stepperObj(y, appObj, resPolObj, jacPolObj);

  using lin_solver_t = solvers::iterative::EigenIterative<
    solvers::linear::iterative::Bicgstab, jac_t>;
  lin_solver_t linSolverObj;
  pressio::solvers::NewtonRaphson<double, lin_solver_t> solverO(linSolverObj);

  // integrate in time
  ::pressio::ode::types::step_t nSteps = 2;
  double dt = 0.01;
  ode::integrateNSteps(stepperObj, y, 0.0, dt, 2, solverO);
  std::cout << std::setprecision(14) << *y.data() << "\n";

  // appObj.analyticAdvanceBackEulerNSteps(dt, nSteps);
  // EXPECT_DOUBLE_EQ(y[0], appObj.y[0]);
  // EXPECT_DOUBLE_EQ(y[1], appObj.y[1]);
  // EXPECT_DOUBLE_EQ(y[2], appObj.y[2]);
}
