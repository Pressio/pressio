
#ifndef ODE_IMPLICIT_METHODS_IMPLICIT_STEPPER_IMPL_HPP_
#define ODE_IMPLICIT_METHODS_IMPLICIT_STEPPER_IMPL_HPP_

#include "ode_implicit_stepper_euler_impl.hpp"
#include "ode_implicit_stepper_bdf2_impl.hpp"
#include "ode_implicit_stepper_arbitrary_impl.hpp"

namespace pressio{ namespace ode{ namespace implicitmethods{ namespace impl{

template<typename tag, typename ...Args>
struct Composer;

////////////////////////////////////////
/// BDF1 stepper
////////////////////////////////////////

// potential cases
// 1. Composer<bdf1, state_t, res_t, jac_t, system_t>;
// 2. Composer<bdf1, state_t, res_t, jac_t, system_t, res_pol, jac_pol>;
// 3. Composer<bdf1, state_t, res_t, jac_t, system_t, void res_pol, jac_pol>;

// 1.
// when we have standard policies, the system must be checked for API
// because the standard policies call some specific methods on the system,
// so the system must meet some requirements
template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename system_type
  >
struct Composer<
  ::pressio::ode::implicitmethods::Euler,
  mpl::enable_if_t<
    ::pressio::ode::concepts::implicit_state<state_type>::value and
    ::pressio::ode::concepts::implicit_residual<residual_type>::value and
    ::pressio::ode::concepts::implicit_jacobian<jacobian_type>::value and
    ::pressio::ode::concepts::continuous_time_implicit_system<system_type>::value
  >,
  state_type, residual_type, jacobian_type, system_type>
{
  using scalar_t = typename ::pressio::containers::details::traits<state_type>::scalar_t;
  static_assert(::pressio::containers::predicates::are_scalar_compatible<state_type,
    residual_type, jacobian_type>::value,
    "state, residual and jacobian are not scalar compatible ");

  using residual_policy_t = ::pressio::ode::implicitmethods::policy::ResidualStandardPolicy<
    state_type, system_type, residual_type>;
  using jacobian_policy_t = ::pressio::ode::implicitmethods::policy::JacobianStandardPolicy<
    state_type, system_type, jacobian_type>;
  using type = StepperBDF1<scalar_t, state_type, residual_type, jacobian_type,
        system_type, residual_policy_t, jacobian_policy_t>;
};

// 2.
// when we have custom policies, we do not check the system, because in general,
// one could not use it since the stepper only knows about the policies, and these can themselves
// contain the system object. Pressio in this case, does not need to checck if system has certain requirements.
// it is user's responsibility to make sure the policies are doing the right thing.
template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename system_type,
  typename residual_policy_type,
  typename jacobian_policy_type
  >
struct Composer<
  ::pressio::ode::implicitmethods::Euler,
  mpl::enable_if_t<
    ::pressio::ode::concepts::implicit_state<state_type>::value and
    ::pressio::ode::concepts::implicit_residual<residual_type>::value and
    ::pressio::ode::concepts::implicit_jacobian<jacobian_type>::value 
    and
    ::pressio::ode::concepts::implicit_euler_residual_policy<
      residual_policy_type, state_type, residual_type, system_type,
      typename ::pressio::containers::details::traits<state_type>::scalar_t>::value 
    and
    ::pressio::ode::concepts::implicit_euler_jacobian_policy<
      jacobian_policy_type, state_type, jacobian_type, system_type,
      typename ::pressio::containers::details::traits<state_type>::scalar_t>::value
  >,
  state_type, residual_type, jacobian_type, system_type, residual_policy_type, jacobian_policy_type>
{
  using scalar_t = typename ::pressio::containers::details::traits<state_type>::scalar_t;
  static_assert(::pressio::containers::predicates::are_scalar_compatible<state_type, residual_type, jacobian_type>::value,
		"state, residual and jacobian are not scalar compatible ");
  using type = StepperBDF1<scalar_t, state_type, residual_type, jacobian_type,
        system_type, residual_policy_type, jacobian_policy_type>;
};

// 3. when we pass a void before policies in place of aux stepper
// when we have custom policies, we do not check the system, because in general,
// one could not use it since the stepper only knows about the policies, and these can themselves
// contain the system object. Pressio in this case, does not need to checck if system has certain requirements.
// it is user's responsibility to make sure the policies are doing the right thing.
template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename system_type,
  typename residual_policy_type,
  typename jacobian_policy_type>
struct Composer<
  ::pressio::ode::implicitmethods::Euler,
  mpl::enable_if_t<
    ::pressio::ode::concepts::implicit_state<state_type>::value and
    ::pressio::ode::concepts::implicit_residual<residual_type>::value and
    ::pressio::ode::concepts::implicit_jacobian<jacobian_type>::value 
    and
    ::pressio::ode::concepts::implicit_euler_residual_policy<
      residual_policy_type, state_type, residual_type, system_type,
      typename ::pressio::containers::details::traits<state_type>::scalar_t>::value 
    and
    ::pressio::ode::concepts::implicit_euler_jacobian_policy<
      jacobian_policy_type, state_type, jacobian_type, system_type,
      typename ::pressio::containers::details::traits<state_type>::scalar_t>::value
  >,
  state_type, residual_type, jacobian_type, system_type, void, residual_policy_type, jacobian_policy_type>
{
  using scalar_t = typename ::pressio::containers::details::traits<state_type>::scalar_t;
  static_assert(::pressio::containers::predicates::are_scalar_compatible<state_type, residual_type, jacobian_type>::value,
		"state, residual and jacobian are not scalar compatible ");
  using type = StepperBDF1<scalar_t, state_type, residual_type, jacobian_type,
        system_type, residual_policy_type, jacobian_policy_type>;
};



// Composer<bdf2, state_t, res_t, jac_t, system_t, aux_stepper_t>;
// Composer<bdf2, state_t, res_t, jac_t, system_t, aux_stepper_t, res_pol, jac_pol>;
template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename system_type,
  typename aux_stepper_t
  >
struct Composer<
  ::pressio::ode::implicitmethods::BDF2,
  mpl::enable_if_t<
    ::pressio::ode::concepts::implicit_state<state_type>::value and
    ::pressio::ode::concepts::implicit_residual<residual_type>::value and
    ::pressio::ode::concepts::implicit_jacobian<jacobian_type>::value and
    ::pressio::ode::concepts::continuous_time_implicit_system<system_type>::value and
    ::pressio::ode::concepts::auxiliary_stepper_for_bdf2<aux_stepper_t>::value
  >,
  state_type, residual_type, jacobian_type, system_type, aux_stepper_t>
{
  using scalar_t = typename ::pressio::containers::details::traits<state_type>::scalar_t;
  static_assert(::pressio::containers::predicates::are_scalar_compatible<state_type,
    residual_type, jacobian_type>::value,
    "state, residual and jacobian are not scalar compatible ");

  using residual_policy_t = ::pressio::ode::implicitmethods::policy::ResidualStandardPolicy<
    state_type, system_type, residual_type>;
  using jacobian_policy_t = ::pressio::ode::implicitmethods::policy::JacobianStandardPolicy<
    state_type, system_type, jacobian_type>;
  using type = StepperBDF2<scalar_t, state_type, residual_type, jacobian_type,
        system_type, aux_stepper_t, residual_policy_t, jacobian_policy_t>;
};

template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename system_type,
  typename aux_stepper_t,
  typename residual_policy_type,
  typename jacobian_policy_type
  >
struct Composer<
  ::pressio::ode::implicitmethods::BDF2,
  mpl::enable_if_t<
    ::pressio::ode::concepts::implicit_state<state_type>::value and
    ::pressio::ode::concepts::implicit_residual<residual_type>::value and
    ::pressio::ode::concepts::implicit_jacobian<jacobian_type>::value and
    ::pressio::ode::concepts::auxiliary_stepper_for_bdf2<aux_stepper_t>::value 
    and
    ::pressio::ode::concepts::implicit_bdf2_residual_policy<
      residual_policy_type, state_type, residual_type, system_type,
      typename ::pressio::containers::details::traits<state_type>::scalar_t>::value 
    and
    ::pressio::ode::concepts::implicit_bdf2_jacobian_policy<
      jacobian_policy_type, state_type, jacobian_type, system_type,
      typename ::pressio::containers::details::traits<state_type>::scalar_t>::value
  >,
  state_type, residual_type, jacobian_type, system_type, aux_stepper_t, residual_policy_type, jacobian_policy_type>
{
  using scalar_t = typename ::pressio::containers::details::traits<state_type>::scalar_t;
  static_assert(::pressio::containers::predicates::are_scalar_compatible<state_type, residual_type, jacobian_type>::value,
		"state, residual and jacobian are not scalar compatible ");

  using type = StepperBDF2<scalar_t, state_type, residual_type, jacobian_type,
			   system_type, aux_stepper_t, residual_policy_type, jacobian_policy_type>;
};



////////////////////////////////////////
/// composer Arbitrary stepper
////////////////////////////////////////

// Composer< state_t, res_t, jac_t, system_t, OrderSetter, TotNumStatesSetter>;
template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename system_type,
  typename order_setter_t,
  typename tot_n_setter_t
  >
struct Composer<
  ::pressio::ode::implicitmethods::Arbitrary,
  mpl::enable_if_t<
    ::pressio::ode::concepts::implicit_state<state_type>::value and
    ::pressio::ode::concepts::implicit_residual<residual_type>::value and
    ::pressio::ode::concepts::implicit_jacobian<jacobian_type>::value and
    ::pressio::ode::concepts::discrete_time_system_implicit_stepping<system_type>::value and
    ::pressio::ode::predicates::IsStepperOrderSetter<order_setter_t>::value and
    ::pressio::ode::predicates::IsStepperTotalNumStatesSetter<tot_n_setter_t>::value
  >,
  state_type, residual_type, jacobian_type, system_type, order_setter_t, tot_n_setter_t>
{
  using scalar_t = typename ::pressio::containers::details::traits<state_type>::scalar_t;
  static_assert(::pressio::containers::predicates::are_scalar_compatible<state_type,
    residual_type, jacobian_type>::value,
    "state, residual and jacobian are not scalar compatible ");

  using residual_policy_t = ::pressio::ode::implicitmethods::policy::ResidualStandardPolicy<
    state_type, system_type, residual_type>;
  using jacobian_policy_t = ::pressio::ode::implicitmethods::policy::JacobianStandardPolicy<
    state_type, system_type, jacobian_type>;

  using type = StepperArbitrary<scalar_t, state_type, residual_type, jacobian_type,
        system_type, order_setter_t, tot_n_setter_t, residual_policy_t, jacobian_policy_t>;
};

// Composer< state_t, res_t, jac_t, system_t, OrderSetter, TotNumStatesSetter, res_pol, jac_pol>;
template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename system_type,
  typename order_setter_t,
  typename tot_n_setter_t,
  typename residual_policy_type,
  typename jacobian_policy_type  
>
struct Composer<
  ::pressio::ode::implicitmethods::Arbitrary,
  mpl::enable_if_t<
    ::pressio::ode::concepts::implicit_state<state_type>::value and
    ::pressio::ode::concepts::implicit_residual<residual_type>::value and
    ::pressio::ode::concepts::implicit_jacobian<jacobian_type>::value and
    ::pressio::ode::predicates::IsStepperOrderSetter<order_setter_t>::value and
    ::pressio::ode::predicates::IsStepperTotalNumStatesSetter<tot_n_setter_t>::value
    and
    ::pressio::ode::concepts::implicit_residual_policy<
      residual_policy_type, ::pressio::ode::implicitmethods::Arbitrary, tot_n_setter_t::value - 1, 
      state_type, residual_type, system_type,
      typename ::pressio::containers::details::traits<state_type>::scalar_t>::value 
    and
    ::pressio::ode::concepts::implicit_jacobian_policy<
      jacobian_policy_type, ::pressio::ode::implicitmethods::Arbitrary, tot_n_setter_t::value - 1, 
      state_type, jacobian_type, system_type,
      typename ::pressio::containers::details::traits<state_type>::scalar_t>::value
  >,
  state_type, residual_type, jacobian_type, system_type, order_setter_t, tot_n_setter_t, 
  residual_policy_type, jacobian_policy_type>
{
  using scalar_t = typename ::pressio::containers::details::traits<state_type>::scalar_t;
  static_assert(::pressio::containers::predicates::are_scalar_compatible<state_type,
    residual_type, jacobian_type>::value,
    "state, residual and jacobian are not scalar compatible ");

  using type = StepperArbitrary<scalar_t, state_type, residual_type, jacobian_type,
        system_type, order_setter_t, tot_n_setter_t, residual_policy_type, jacobian_policy_type>;
};


}}}} // end namespace pressio::ode::explicitmethods::impl
#endif
