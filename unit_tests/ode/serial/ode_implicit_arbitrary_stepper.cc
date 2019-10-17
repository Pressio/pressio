
#include <gtest/gtest.h>
#include "CONTAINERS_ALL"
#include "ODE_ALL"
#include "reference_apps_for_testing.hpp"

template<typename state_type, typename system_type, typename residual_type>
class ResidualPolicy
  : public ::pressio::ode::policy::ImplicitResidualPolicyBase<
  ResidualPolicy<state_type, system_type, residual_type>
  >
{
public:

  static constexpr auto stepper_order = 3;
  static constexpr auto num_aux_states = 1;

  void operator()(const state_type & y,
		  residual_type & R,
		  const std::array<state_type, num_aux_states> & oldYs,
		  const system_type & model,
		  double t,
		  double dt,
		  ::pressio::ode::types::step_t step) const{
    // here I would need to compute the time discrete residual
  }

  residual_type operator()(const state_type & y,
  			   const std::array<state_type, num_aux_states> & oldYs,
  			   const system_type & model,
  			   double t,
  			   double dt,
			   ::pressio::ode::types::step_t step) const{
    // here I would need to compute the time discrete residual
    residual_type R;
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
  static constexpr auto stepper_order = 3;
  static constexpr auto num_aux_states = 1;

  void operator()(const state_type & y,
		  jacobian_type & J,
		  const system_type & model,
		  double t,
		  double dt,
		  ::pressio::ode::types::step_t step) const{
    // here I would need to compute the time discrete version
  }

  jacobian_type operator()(const state_type & y,
  			   const system_type & model,
  			   double t,
  			   double dt,
			   ::pressio::ode::types::step_t step) const{
    // here I would need to compute the time discrete version
    jacobian_type J;
    return J;
  }

};//end class


TEST(ode_implicit, validArbitraryStepperPolicies){
  using namespace pressio;

  using app_t = ode::testing::refAppForImpEigen;
  using nstate_t = typename app_t::state_type;
  using nvel_t = typename app_t::velocity_type;
  using njac_t = typename app_t::jacobian_type;
  using state_t = containers::Vector<nstate_t>;
  using res_t = containers::Vector<nvel_t>;
  using jac_t = containers::Matrix<njac_t>;

  using residual_policy_t = ResidualPolicy<state_t, app_t, res_t>;
  static_assert
    (ode::meta::is_legitimate_residual_policy_for_implicit_arbitrary_stepper<
     residual_policy_t, state_t, res_t, app_t, double>::value, "");

  using jacobian_policy_t = JacobianPolicy<state_t, app_t, jac_t>;
  static_assert
    (ode::meta::is_legitimate_jacobian_policy_for_implicit_arbitrary_stepper<
     jacobian_policy_t, state_t, jac_t, app_t, double>::value, "");
}


TEST(ode_implicit, validArbitraryStepper){
  using namespace pressio;

  using app_t = ode::testing::refAppForImpEigen;
  using nstate_t = typename app_t::state_type;
  using nvel_t = typename app_t::velocity_type;
  using njac_t = typename app_t::jacobian_type;
  using state_t = containers::Vector<nstate_t>;
  using res_t = containers::Vector<nvel_t>;
  using jac_t = containers::Matrix<njac_t>;

  using residual_policy_t = ResidualPolicy<state_t, app_t, res_t>;
  static_assert
    (ode::meta::is_legitimate_residual_policy_for_implicit_arbitrary_stepper<
     residual_policy_t, state_t, res_t, app_t, double>::value, "");

  using jacobian_policy_t = JacobianPolicy<state_t, app_t, jac_t>;
  static_assert
    (ode::meta::is_legitimate_jacobian_policy_for_implicit_arbitrary_stepper<
     jacobian_policy_t, state_t, jac_t, app_t, double>::value, "");

  using stepper_t = ode::ImplicitStepper<
    ode::ImplicitEnum::Arbitrary,
    state_t, res_t, jac_t, app_t,
    residual_policy_t, jacobian_policy_t>;

  using traits = ode::details::traits<stepper_t>;
  static_assert( traits::is_implicit, "");

  static_assert( traits::order_value == 3, "");
  static_assert( traits::numAuxStates == 1, "");

  ::testing::StaticAssertTypeEq<typename
  				traits::state_t, state_t>();
  ::testing::StaticAssertTypeEq<typename
  				traits::residual_t,res_t>();
  ::testing::StaticAssertTypeEq<typename
  				traits::jacobian_t,jac_t>();
  ::testing::StaticAssertTypeEq<typename
  				traits::scalar_t,double>();
  ::testing::StaticAssertTypeEq<typename
  				traits::system_t,app_t>();
}
