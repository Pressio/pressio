
#include <gtest/gtest.h>
#include "pressio_ode.hpp"
#include "../reference_apps_for_testing.hpp"

template<typename state_type, typename system_type, typename residual_type>
class ResidualPolicy{

public:
  template <typename prev_states_type>
  void compute(const state_type & y,
		  const prev_states_type & oldYs,
		  const system_type & model,
		  const double & t,
		  const double & dt,
		  ::pressio::ode::types::step_t step,
		  residual_type & R,
      ::pressio::Norm normKind,
      double & normValue) const
  {
    // here I would need to compute the time discrete residual
  }

  residual_type create(const system_type & model) const{
    return residual_type();
  }

};//end class


template<typename state_type, typename system_type, typename jacobian_type>
class JacobianPolicy{

public:
  template <typename prev_states_type>
  void compute(const state_type & y,
		  const prev_states_type & oldYs,
		  const system_type & model,
		  const double &  t,
		  const double &  dt,
		  ::pressio::ode::types::step_t step,
		  jacobian_type & J) const{
    // here I would need to compute the time discrete version
  }

  jacobian_type create(const system_type & model) const{
    return jacobian_type();
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
    (ode::meta::admissible_residual_policy_for_implicit_arbitrary_stepper<
     residual_policy_t, 1, state_t, res_t, app_t, double>::value, "");

  using jacobian_policy_t = JacobianPolicy<state_t, app_t, jac_t>;
  static_assert
    (ode::meta::admissible_jacobian_policy_for_implicit_arbitrary_stepper<
     jacobian_policy_t, 1, state_t, jac_t, app_t, double>::value, "");
}


// TEST(ode_implicit, validArbitraryStepper){
//   using namespace pressio;

//   using app_t	 = ode::testing::refAppForArbitraryImpl;
//   using nstate_t = typename app_t::state_type;
//   using nres_t   = typename app_t::residual_type;
//   using njac_t	 = typename app_t::jacobian_type;
//   using state_t  = containers::Vector<nstate_t>;
//   using res_t	 = containers::Vector<nres_t>;
//   using jac_t    = containers::Matrix<njac_t>;

//   using stepper_order = ode::types::StepperOrder<1>;
//   using stepper_n_states = ode::types::StepperTotalNumberOfStates<2>;

//   using stepper_t = ode::ImplicitStepper<
//     ode::implicitmethods::Arbitrary,
//     state_t, res_t, jac_t, app_t,
//     stepper_order, stepper_n_states>;

//   using traits = ode::details::traits<stepper_t>;
//   static_assert( traits::is_implicit, "");

//   static_assert( traits::order_value == 1, "");
//   // numAuxStates = 1 because it is one less than the total stepper states
//   static_assert( traits::numAuxStates == 1, "");

//   ::testing::StaticAssertTypeEq<typename
//   				traits::state_t, state_t>();
//   ::testing::StaticAssertTypeEq<typename
//   				traits::residual_t,res_t>();
//   ::testing::StaticAssertTypeEq<typename
//   				traits::jacobian_t,jac_t>();
//   ::testing::StaticAssertTypeEq<typename
//   				traits::scalar_t,double>();
//   ::testing::StaticAssertTypeEq<typename
//   				traits::system_t,app_t>();
// }
