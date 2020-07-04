
#ifndef ODE_IMPLICIT_METHODS_IMPLICIT_STEPPER_IMPL_HPP_
#define ODE_IMPLICIT_METHODS_IMPLICIT_STEPPER_IMPL_HPP_

#include "ode_implicit_stepper_euler_impl.hpp"
#include "ode_implicit_stepper_bdf2_impl.hpp"
#include "ode_implicit_stepper_arbitrary_impl.hpp"

namespace pressio{ namespace ode{ namespace implicitmethods{ namespace impl{

template<typename tag, typename ...Args> 
struct RegStepperOptionals;

// if we fall here, no ud_ops and standard policies
template<
  typename tag, typename state_t, 
  typename residual_t, typename jacobian_t, 
  typename system_t, typename scalar_t
  >
struct RegStepperOptionals<
  tag,
  mpl::enable_if_t<
   std::is_same<tag, ::pressio::ode::implicitmethods::Euler>::value or 
   std::is_same<tag, ::pressio::ode::implicitmethods::BDF2>::value    
  >,   
  state_t, residual_t, jacobian_t, system_t, scalar_t
>
{
  using residual_policy_type = ::pressio::ode::implicitmethods::policy::ResidualStandardPolicy<
    state_t, system_t, residual_t>;
  using jacobian_policy_type = ::pressio::ode::implicitmethods::policy::JacobianStandardPolicy<
    state_t, system_t, jacobian_t>;
  using ud_ops_t = void;
};

// if we fall here, no ud_ops and but custom policies
template<
  typename state_t, typename residual_t, typename jacobian_t, typename system_t, 
  typename scalar_t, typename res_pol, typename jac_pol
  >
struct RegStepperOptionals<
 ::pressio::ode::implicitmethods::Euler,
  mpl::enable_if_t<
   ::pressio::ode::meta::admissible_implicit_euler_residual_policy_regular_stepper<
    res_pol, state_t, residual_t, system_t, scalar_t>::value 
    and 
   ::pressio::ode::meta::admissible_implicit_euler_jacobian_policy_regular_stepper<
    jac_pol, state_t, jacobian_t, system_t, scalar_t>::value
  >,   
  state_t, residual_t, jacobian_t, system_t, scalar_t, res_pol, jac_pol
>
{
  using residual_policy_type = res_pol;
  using jacobian_policy_type = jac_pol;
  using ud_ops_t = void;
};

// if we fall here, no ud_ops and but custom policies
template<
  typename state_t, typename residual_t, typename jacobian_t, typename system_t, 
  typename scalar_t, typename res_pol, typename jac_pol
  >
struct RegStepperOptionals<
 ::pressio::ode::implicitmethods::BDF2,
  mpl::enable_if_t<
   ::pressio::ode::meta::admissible_implicit_bdf2_residual_policy_regular_stepper<
    res_pol, state_t, residual_t, system_t, scalar_t>::value 
    and 
   ::pressio::ode::meta::admissible_implicit_bdf2_jacobian_policy_regular_stepper<
    jac_pol, state_t, jacobian_t, system_t, scalar_t>::value
  >,   
  state_t, residual_t, jacobian_t, system_t, scalar_t, res_pol, jac_pol
>
{
  using residual_policy_type = res_pol;
  using jacobian_policy_type = jac_pol;
  using ud_ops_t = void;
};



////////////////////////////////////////
/// composer 
////////////////////////////////////////

template<typename tag, typename ...Args> 
struct Composer;


////////////////////////////////////////
/// composer regular stepper
////////////////////////////////////////

// Composer<bdf1, state_t, res_t, jac_t, system_t>;
// Composer<bdf1, state_t, res_t, jac_t, system_t, res_pol, jac_pol>;
template<
  typename state_type, 
  typename residual_type, 
  typename jacobian_type, 
  typename system_type, 
  typename ...Args>
struct Composer<
  ::pressio::ode::implicitmethods::Euler, 
  mpl::enable_if_t<
    ::pressio::ode::meta::legitimate_implicit_state_type<state_type>::value and
    ::pressio::ode::meta::legitimate_implicit_residual_type<residual_type>::value and 
    ::pressio::ode::meta::legitimate_jacobian_type<jacobian_type>::value and     
    ::pressio::ode::meta::admissible_system_implicit_ode_regular_stepper<system_type>::value
  >, 
  state_type, residual_type, jacobian_type, system_type, Args...>
{ 
  using scalar_t = typename ::pressio::containers::details::traits<state_type>::scalar_t;
  static_assert(::pressio::containers::meta::are_scalar_compatible<state_type, 
    residual_type, jacobian_type>::value, 
    "state, residual and jacobian are not scalar compatible ");

  using prop = RegStepperOptionals<
     ::pressio::ode::implicitmethods::Euler, 
     void, state_type, residual_type, jacobian_type, system_type, scalar_t, Args...>;
  using residual_policy_t = typename prop::residual_policy_type;
  using jacobian_policy_t = typename prop::jacobian_policy_type;

  using type = StepperBDF1<scalar_t, state_type, residual_type, jacobian_type, 
        system_type, residual_policy_t, jacobian_policy_t>;   
};


// Composer<bdf2, state_t, res_t, jac_t, system_t, aux_stepper_t>;
// Composer<bdf2, state_t, res_t, jac_t, system_t, aux_stepper_t, res_pol, jac_pol>;
template<
  typename state_type, 
  typename residual_type, 
  typename jacobian_type, 
  typename system_type, 
  typename aux_stepper_t, 
  typename ...Args>
struct Composer<
  ::pressio::ode::implicitmethods::BDF2, 
  mpl::enable_if_t<
    ::pressio::ode::meta::legitimate_implicit_state_type<state_type>::value and
    ::pressio::ode::meta::legitimate_implicit_residual_type<residual_type>::value and 
    ::pressio::ode::meta::legitimate_jacobian_type<jacobian_type>::value and     
    ::pressio::ode::meta::admissible_system_implicit_ode_regular_stepper<system_type>::value 
    // ::pressio::ode::meta::legitimate_auxiliary_stepper_for_implicit_ode<aux_stepper_t>::value
  >, 
  state_type, residual_type, jacobian_type, system_type, aux_stepper_t, Args...>
{ 
  using scalar_t = typename ::pressio::containers::details::traits<state_type>::scalar_t;
  static_assert(::pressio::containers::meta::are_scalar_compatible<state_type, 
    residual_type, jacobian_type>::value, 
    "state, residual and jacobian are not scalar compatible ");

  using prop = RegStepperOptionals<
     ::pressio::ode::implicitmethods::Euler, 
     void, state_type, residual_type, jacobian_type, system_type, scalar_t, Args...>;
  using residual_policy_t = typename prop::residual_policy_type;
  using jacobian_policy_t = typename prop::jacobian_policy_type;

  using type = StepperBDF2<scalar_t, state_type, residual_type, jacobian_type, 
        system_type, aux_stepper_t, residual_policy_t, jacobian_policy_t>;  
};


////////////////////////////////////////
/// composer Arbitrary stepper
////////////////////////////////////////


// template<typename res_pol, typename jac_pol>
// struct ArbStepperPropUnpackArgs<
//   mpl::enable_if_t<
//   admissible_residual_policy_for_implicit_arbitrary_stepper<res_pol>::value and 
//   admissible_jacobian_policy_for_implicit_arbitrary_stepper<jac_pol>::value 
//   >, res_pol, jac_pol
//   >
// {
//   using residual_policy_t = void;
//   using jacobian_policy_t = void;
// };

template<typename state_t, typename ...Args> 
struct ArbStepperOptionals;

// if we fall here, no ud_ops and standard policies
template<
  typename state_t, typename residual_t, typename jacobian_t, 
  typename system_t, typename scalar_t
  >
struct ArbStepperOptionals<
  state_t, residual_t, jacobian_t, system_t, scalar_t
>
{
  using residual_policy_type = ::pressio::ode::implicitmethods::policy::ResidualStandardPolicyForArbitraryStepper<
    state_t, system_t, residual_t>;
  using jacobian_policy_type = ::pressio::ode::implicitmethods::policy::JacobianStandardPolicyForArbitraryStepper<
    state_t, system_t, jacobian_t>;
  using ud_ops_t = void;
};

// Composer< state_t, res_t, jac_t, system_t, OrderSetter, TotNumStatesSetter>;
// Composer< state_t, res_t, jac_t, system_t, OrderSetter, TotNumStatesSetter, res_pol, jac_pol>;
template<
  typename state_type, 
  typename residual_type,  
  typename jacobian_type, 
  typename system_type, 
  typename order_setter_t, 
  typename tot_n_setter_t, 
  typename ...Args
  >
struct Composer<
  ::pressio::ode::implicitmethods::Arbitrary, 
  mpl::enable_if_t<
    ::pressio::ode::meta::legitimate_implicit_state_type<state_type>::value and
    ::pressio::ode::meta::legitimate_implicit_residual_type<residual_type>::value and 
    ::pressio::ode::meta::legitimate_jacobian_type<jacobian_type>::value and
    ::pressio::ode::meta::admissible_system_implicit_ode_arbitrary_stepper<system_type>::value and
    ::pressio::ode::meta::IsStepperOrderSetter<order_setter_t>::value and   
    ::pressio::ode::meta::IsStepperTotalNumStatesSetter<tot_n_setter_t>::value 
  >, 
  state_type, residual_type, jacobian_type, system_type, order_setter_t, tot_n_setter_t, Args...>
{ 
  using scalar_t = typename ::pressio::containers::details::traits<state_type>::scalar_t;
  static_assert(::pressio::containers::meta::are_scalar_compatible<state_type, 
    residual_type, jacobian_type>::value, 
    "state, residual and jacobian are not scalar compatible ");

  using prop = ArbStepperOptionals<state_type, residual_type, jacobian_type, system_type, scalar_t, Args...>;
  using residual_policy_t = typename prop::residual_policy_type;
  using jacobian_policy_t = typename prop::jacobian_policy_type;

  using type = StepperArbitrary<scalar_t, state_type, residual_type, jacobian_type, 
        system_type, order_setter_t, tot_n_setter_t, residual_policy_t, jacobian_policy_t>;
};


}}}} // end namespace pressio::ode::explicitmethods::impl
#endif
