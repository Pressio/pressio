
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_TRAITS_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_TRAITS_HPP_

#include "../../ode_fwd.hpp"
#include "../policies/meta/ode_find_if_legitimate_implicit_residual_policy.hpp"
#include "../policies/meta/ode_find_if_legitimate_implicit_jacobian_policy.hpp"
#include "../../meta/ode_is_valid_user_defined_ops_for_implicit_ode.hpp"

namespace rompp{ namespace ode{ namespace details{

template <typename T, typename = void>
struct ScalarHelper{
  static constexpr bool value = false;
  using type = void;
};

template <typename T>
struct ScalarHelper<
  T,
  mpl::enable_if_t<
    ::rompp::algebra::meta::is_wrapper<T>::value and
    ::rompp::algebra::details::traits<T>::wrapped_package_identifier
    != ::rompp::algebra::details::WrappedPackageIdentifier::Arbitrary
    >
  >{
  static constexpr bool value = true;
  using type = typename ::rompp::algebra::details::traits<T>::scalar_t;
};

template <typename T>
struct ScalarHelper<
  T,
  mpl::enable_if_t<
    ::rompp::algebra::meta::is_wrapper<T>::value and
    ::rompp::algebra::details::traits<T>::wrapped_package_identifier
    == ::rompp::algebra::details::WrappedPackageIdentifier::Arbitrary
    >
  >{
  static constexpr bool value = false;
  using type = void;
};
//-------------------------------------------------------------------


template <
  typename model_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t,
  typename enable = void
  >
struct StdPoliciesPicker;

template <
  typename model_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t
  >
struct StdPoliciesPicker<
  model_t, state_t, residual_t, jacobian_t
#ifdef HAVE_PYBIND11  
  , mpl::enable_if_t<
    mpl::not_same<model_t, pybind11::object>::value
    >
#endif    
  >
{
  using standard_res_policy_t = policy::ImplicitResidualStandardPolicy<
    state_t, model_t, residual_t>;
  using standard_jac_policy_t = policy::ImplicitJacobianStandardPolicy<
    state_t, model_t, jacobian_t>;
};

#ifdef HAVE_PYBIND11  
template <
  typename model_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t
  >
struct StdPoliciesPicker<
  model_t, state_t, residual_t, jacobian_t,
  mpl::enable_if_t<
    mpl::is_same<model_t, pybind11::object>::value
    >
  >
{
  using standard_res_policy_t = policy::ImplicitResidualStandardPolicyPybind11<
    state_t, model_t, residual_t>;
  using standard_jac_policy_t = policy::ImplicitJacobianStandardPolicyPybind11<
    state_t, model_t, jacobian_t>;
};
#endif
//------------------------------------------------------------------


template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename model_type,
  typename ...Args
  >
struct traits<
  ImplicitStepper<
    ImplicitEnum::Euler,
    state_type, residual_type,
    jacobian_type, model_type,
    Args...>
  > {

  using stepper_t =   ImplicitStepper< ImplicitEnum::Euler,
				       state_type, residual_type,
				       jacobian_type, model_type,
				       Args...>;
  using this_t = traits<stepper_t>;

  static constexpr ode::ImplicitEnum enum_id = ode::ImplicitEnum::Euler;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;

  using state_t		  = state_type;
  using residual_t	  = residual_type;
  using jacobian_t	  = jacobian_type;
  using model_t		  = model_type;
  using aux_stepper_t	  = void;

  static constexpr unsigned int order_value = 1;
  static constexpr unsigned int steps = 1;

  // check if scalar is provided in Args
  using ic0 = ::rompp::mpl::variadic::find_if_unary_pred_t<std::is_floating_point, Args...>;
  using scalar_from_args = ::rompp::mpl::variadic::at_or_t<void, ic0::value, Args...>;
  // check if state is a algebra wrapper, and if so get its scalar_type
  using scalar_type_from_traits = typename ScalarHelper<state_type>::type;
  // decide which to pick
  using scalar_t = typename std::conditional<
    std::is_void<scalar_from_args>::value,
    scalar_type_from_traits, scalar_from_args>::type;

  static_assert( std::is_floating_point<scalar_t>::value,
  		 "I cannot guess the scalar_type because it is not found in templates and the state_type used is not a algebra wrapper. If you are using custom data structures that do not have wrappers in the algebra, pass scalar as a template.");


  // standard policies (only used if not passed a user-defined policy)
  using policy_picker = StdPoliciesPicker<model_t, state_t, residual_t, jacobian_t>;
  using standard_res_policy_t = typename policy_picker::standard_res_policy_t;
  using standard_jac_policy_t = typename policy_picker::standard_jac_policy_t;

  // check Args for a user-defined admissible residual policy
  using ic1 = ::rompp::ode::meta::find_if_legitimate_implicit_residual_policy_t<
    this_t::enum_id, this_t::steps, state_t, residual_t, model_t, scalar_t,
    Args...>;
  using residual_policy_t = ::rompp::mpl::variadic::at_or_t
    <standard_res_policy_t, ic1::value, Args...>;

  // check Args for a user-defined admissible jacobian policy
  using ic2 = ::rompp::ode::meta::find_if_legitimate_implicit_jacobian_policy_t<
    this_t::enum_id, state_t, jacobian_t, model_t, scalar_t,
    Args...>;
  using jacobian_policy_t = ::rompp::mpl::variadic::at_or_t
    <standard_jac_policy_t, ic2::value, Args...>;
};


template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename model_type,
  typename ...Args
  >
struct traits<
  ImplicitStepper<
    ImplicitEnum::BDF2,
    state_type, residual_type,
    jacobian_type, model_type,
    Args...>
  > {

  using stepper_t =   ImplicitStepper< ImplicitEnum::BDF2,
				       state_type, residual_type,
				       jacobian_type, model_type,
				       Args...>;
  using this_t = traits<stepper_t>;

  static constexpr ode::ImplicitEnum enum_id = ode::ImplicitEnum::BDF2;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;

  using state_t		  = state_type;
  using residual_t	  = residual_type;
  using jacobian_t	  = jacobian_type;
  using model_t		  = model_type;

  static constexpr unsigned int order_value = 2;
  static constexpr unsigned int steps = 2;

  // check if scalar is provided in Args
  using ic0 = ::rompp::mpl::variadic::find_if_unary_pred_t<std::is_floating_point, Args...>;
  using scalar_from_args = ::rompp::mpl::variadic::at_or_t<void, ic0::value, Args...>;
  // check if state is a algebra wrapper, and if so get its scalar_type
  using scalar_type_from_traits = typename ScalarHelper<state_type>::type;
  // decide which to pick
  using scalar_t = typename std::conditional<
    std::is_void<scalar_from_args>::value,
    scalar_type_from_traits, scalar_from_args>::type;

  static_assert( std::is_floating_point<scalar_t>::value,
  		 "I cannot guess the scalar_type because it is not found in templates and the state_type used is not a algebra wrapper. If you are using custom data structures that do not have wrappers in the algebra, pass scalar as a template.");


  // for BDF2 the user has to pass an auxiliary stepper
  using ic1 = ::rompp::mpl::variadic::find_if_binary_pred_t<
    stepper_t, ::rompp::ode::meta::is_legitimate_auxiliary_stepper, Args...>;
  using aux_stepper_t = ::rompp::mpl::variadic::at_or_t<void, ic1::value, Args...>;

  // // standard policies (only used if user-defined policies not passed)
  using policy_picker = StdPoliciesPicker<model_t, state_t, residual_t, jacobian_t>;
  using standard_res_policy_t = typename policy_picker::standard_res_policy_t;
  using standard_jac_policy_t = typename policy_picker::standard_jac_policy_t;

  // check Args if a user-defined admissible residual policy is passed
  using ic2 = ::rompp::ode::meta::find_if_legitimate_implicit_residual_policy_t<
    this_t::enum_id, this_t::steps, state_t, residual_t, model_t, scalar_t,
    Args...>;
  using residual_policy_t = ::rompp::mpl::variadic::at_or_t
    <standard_res_policy_t, ic2::value, Args...>;

  // check Args if a user-defined admissible jacobian policy is passed
  using ic3 = ::rompp::ode::meta::find_if_legitimate_implicit_jacobian_policy_t<
    this_t::enum_id, state_t, jacobian_t, model_t, scalar_t,
    Args...>;
  using jacobian_policy_t = ::rompp::mpl::variadic::at_or_t
    <standard_jac_policy_t, ic3::value, Args...>;
};

}}}//end namespace rompp::ode::details
#endif
