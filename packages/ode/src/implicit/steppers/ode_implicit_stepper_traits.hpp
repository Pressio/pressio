/*
//@HEADER
// ************************************************************************
//
// ode_implicit_stepper_traits.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_TRAITS_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_TRAITS_HPP_

#include "../../ode_fwd.hpp"
#include "../policies/meta/ode_find_if_legitimate_implicit_residual_policy.hpp"
#include "../policies/meta/ode_find_if_legitimate_implicit_jacobian_policy.hpp"
#include "../policies/meta/ode_find_if_legitimate_residual_policy_for_implicit_arbitrary_stepper.hpp"
#include "../policies/meta/ode_find_if_legitimate_jacobian_policy_for_implicit_arbitrary_stepper.hpp"
#include "../../meta/ode_is_valid_user_defined_ops_for_implicit_ode.hpp"

namespace pressio{ namespace ode{ namespace details{

namespace impl{

template <typename T, typename = void>
struct ScalarHelper{
  static constexpr bool value = false;
  using type = void;
};

template <typename T>
struct ScalarHelper<
  T,
  mpl::enable_if_t<
    ::pressio::containers::meta::is_wrapper<T>::value and
    ::pressio::containers::details::traits<T>::wrapped_package_identifier
    != ::pressio::containers::details::WrappedPackageIdentifier::Arbitrary
    >
  >{
  static constexpr bool value = true;
  using type = typename ::pressio::containers::details::traits<T>::scalar_t;
};

template <typename T>
struct ScalarHelper<
  T,
  mpl::enable_if_t<
    ::pressio::containers::meta::is_wrapper<T>::value and
    ::pressio::containers::details::traits<T>::wrapped_package_identifier
    == ::pressio::containers::details::WrappedPackageIdentifier::Arbitrary
    >
  >{
  static constexpr bool value = false;
  using type = void;
};
//-------------------------------------------------------------------


template <
  typename system_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t,
  typename enable = void
  >
struct StdPoliciesPicker;

template <
  typename system_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t
  >
struct StdPoliciesPicker<
  system_t, state_t, residual_t, jacobian_t
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  , mpl::enable_if_t<
    mpl::not_same<system_t, pybind11::object>::value
    >
#endif
  >
{
  using standard_res_policy_t = policy::ImplicitResidualStandardPolicy<
    state_t, system_t, residual_t>;
  using standard_jac_policy_t = policy::ImplicitJacobianStandardPolicy<
    state_t, system_t, jacobian_t>;
};

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <
  typename system_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t
  >
struct StdPoliciesPicker<
  system_t, state_t, residual_t, jacobian_t,
  mpl::enable_if_t<
    mpl::is_same<system_t, pybind11::object>::value
    >
  >
{
  using standard_res_policy_t = policy::ImplicitResidualStandardPolicyPybind11<
    state_t, system_t, residual_t>;
  using standard_jac_policy_t = policy::ImplicitJacobianStandardPolicyPybind11<
    state_t, system_t, jacobian_t>;
};
#endif

} // end namepsace pressio::ode::details::impl


template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename system_type,
  typename ...Args
  >
struct traits<
  ImplicitStepper<
    ImplicitEnum::Euler,
    state_type, residual_type,
    jacobian_type, system_type,
    Args...>
  > {

  using stepper_t =   ImplicitStepper< ImplicitEnum::Euler,
				       state_type, residual_type,
				       jacobian_type, system_type,
				       Args...>;
  using this_t = traits<stepper_t>;

  static constexpr ode::ImplicitEnum enum_id = ode::ImplicitEnum::Euler;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;

  using state_t		  = state_type;
  using residual_t	  = residual_type;
  using jacobian_t	  = jacobian_type;
  using system_t	  = system_type;
  using aux_stepper_t	  = void;

  static constexpr types::stepper_order_t order_value = 1;
  static constexpr types::stepper_n_states_t numAuxStates = 1;

  // check if scalar is provided in Args
  using ic0 = ::pressio::mpl::variadic::find_if_unary_pred_t<std::is_floating_point, Args...>;
  using scalar_from_args = ::pressio::mpl::variadic::at_or_t<void, ic0::value, Args...>;
  // check if state is a containers wrapper, and if so get its scalar_type
  using scalar_type_from_traits = typename impl::ScalarHelper<state_type>::type;
  // decide which to pick
  using scalar_t = typename std::conditional<
    std::is_void<scalar_from_args>::value,
    scalar_type_from_traits, scalar_from_args>::type;

  static_assert( std::is_floating_point<scalar_t>::value,
  		 "I cannot determine the scalar_type because it is not found \
in the templates args and the state_type used is not a containers wrapper. \
If you are using custom data structures that do not have wrappers in \
the containers, pass scalar as a template.");

  // standard policies (only used if not passed a user-defined policy)
  using policy_picker = impl::StdPoliciesPicker<system_t, state_t, residual_t, jacobian_t>;
  using standard_res_policy_t = typename policy_picker::standard_res_policy_t;
  using standard_jac_policy_t = typename policy_picker::standard_jac_policy_t;

  // check Args for a user-defined admissible residual policy
  using ic1 = ::pressio::ode::meta::find_if_legitimate_implicit_residual_policy_t<
    this_t::enum_id, this_t::numAuxStates, state_t, residual_t, system_t, scalar_t,
    Args...>;
  using residual_policy_t = ::pressio::mpl::variadic::at_or_t
    <standard_res_policy_t, ic1::value, Args...>;

  // check Args for a user-defined admissible jacobian policy
  using ic2 = ::pressio::ode::meta::find_if_legitimate_implicit_jacobian_policy_t<
    this_t::enum_id, state_t, jacobian_t, system_t, scalar_t,
    Args...>;
  using jacobian_policy_t = ::pressio::mpl::variadic::at_or_t
    <standard_jac_policy_t, ic2::value, Args...>;
};


template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename system_type,
  typename ...Args
  >
struct traits<
  ImplicitStepper<
    ImplicitEnum::BDF2,
    state_type, residual_type,
    jacobian_type, system_type,
    Args...>
  > {

  using stepper_t =   ImplicitStepper< ImplicitEnum::BDF2,
				       state_type, residual_type,
				       jacobian_type, system_type,
				       Args...>;
  using this_t = traits<stepper_t>;

  static constexpr ode::ImplicitEnum enum_id = ode::ImplicitEnum::BDF2;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;

  using state_t	   = state_type;
  using residual_t = residual_type;
  using jacobian_t = jacobian_type;
  using system_t   = system_type;

  static constexpr types::stepper_order_t order_value = 2;
  static constexpr types::stepper_n_states_t numAuxStates = 2;

  // check if scalar is provided in Args
  using ic0 = ::pressio::mpl::variadic::find_if_unary_pred_t<std::is_floating_point, Args...>;
  using scalar_from_args = ::pressio::mpl::variadic::at_or_t<void, ic0::value, Args...>;
  // check if state is a containers wrapper, and if so get its scalar_type
  using scalar_type_from_traits = typename impl::ScalarHelper<state_type>::type;
  // decide which to pick
  using scalar_t = typename std::conditional<
    std::is_void<scalar_from_args>::value,
    scalar_type_from_traits, scalar_from_args>::type;

  static_assert( std::is_floating_point<scalar_t>::value,
  		 "I cannot determine the scalar_type because it is not found \
in the templates args and the state_type used is not a containers wrapper. \
If you are using custom data structures that do not have wrappers in \
the containers, pass scalar as a template.");

  // for BDF2 the user has to pass an auxiliary stepper
  using ic1 = ::pressio::mpl::variadic::find_if_binary_pred_t<
    stepper_t, ::pressio::ode::meta::is_legitimate_auxiliary_stepper, Args...>;
  using aux_stepper_t = ::pressio::mpl::variadic::at_or_t<void, ic1::value, Args...>;

  // standard policies (only used if user-defined policies not passed)
  using policy_picker = impl::StdPoliciesPicker<system_t, state_t, residual_t, jacobian_t>;
  using standard_res_policy_t = typename policy_picker::standard_res_policy_t;
  using standard_jac_policy_t = typename policy_picker::standard_jac_policy_t;

  // check Args if a user-defined admissible residual policy is passed
  using ic2 = ::pressio::ode::meta::find_if_legitimate_implicit_residual_policy_t<
    this_t::enum_id, this_t::numAuxStates, state_t, residual_t, system_t, scalar_t,
    Args...>;
  using residual_policy_t = ::pressio::mpl::variadic::at_or_t
    <standard_res_policy_t, ic2::value, Args...>;

  // check Args if a user-defined admissible jacobian policy is passed
  using ic3 = ::pressio::ode::meta::find_if_legitimate_implicit_jacobian_policy_t<
    this_t::enum_id, state_t, jacobian_t, system_t, scalar_t,
    Args...>;
  using jacobian_policy_t = ::pressio::mpl::variadic::at_or_t
    <standard_jac_policy_t, ic3::value, Args...>;
};




template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename system_type,
  typename ...Args
  >
struct traits<
  ImplicitStepper<
    ImplicitEnum::Arbitrary,
    state_type, residual_type,
    jacobian_type, system_type, Args...>
  > {

  using stepper_t =   ImplicitStepper< ImplicitEnum::Arbitrary,
				       state_type, residual_type,
				       jacobian_type, system_type, Args...>;
  using this_t = traits<stepper_t>;

  static constexpr ode::ImplicitEnum enum_id = ode::ImplicitEnum::Arbitrary;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;

  using state_t		  = state_type;
  using residual_t	  = residual_type;
  using jacobian_t	  = jacobian_type;
  using system_t	  = system_type;
  using aux_stepper_t	  = void;

  // standard policies (only used if not passed a user-defined policy)
  using policy_picker = impl::StdPoliciesPicker<system_t, state_t, residual_t, jacobian_t>;
  using standard_res_policy_t = typename policy_picker::standard_res_policy_t;
  using standard_jac_policy_t = typename policy_picker::standard_jac_policy_t;

  //-----------------------------------------------------------
  // check if scalar is provided in Args
  //-----------------------------------------------------------
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<std::is_floating_point, Args...>;
  using scalar_from_args = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  // check if state is a containers wrapper, and if so get its scalar_type
  using scalar_type_from_traits = typename impl::ScalarHelper<state_type>::type;
  // decide which to pick
  using scalar_t = typename std::conditional<std::is_void<scalar_from_args>::value,
					     scalar_type_from_traits,
					     scalar_from_args>::type;

  static_assert( std::is_floating_point<scalar_t>::value,
  		 "I cannot determine the scalar_type because it is not found \
in the templates args and the state_type used is not a containers' wrapper. \
If you are using custom data structures that are not wrappeed in containers types, \
pass a valid scalar type as a template argument to ImplicitStepper.");

  //-----------------------------------------------------------
  // check Args for a user-defined admissible residual policy
  //-----------------------------------------------------------
  // using ic3 = ::pressio::mpl::variadic::find_if_quinary_pred_t<
  //   state_t, residual_t, system_t, scalar_t,
  //   meta::is_legitimate_residual_policy_for_implicit_arbitrary_stepper, Args...>;
  using ic3 = ::pressio::ode::meta::find_if_legitimate_residual_policy_for_implicit_arbitrary_stepper_t<
    state_t, residual_t, system_t, scalar_t, Args...>;
  using residual_policy_t = ::pressio::mpl::variadic::at_or_t<void, ic3::value, Args...>;

  // check Args for a user-defined admissible jacobian policy
  using ic4 = ::pressio::ode::meta::find_if_legitimate_jacobian_policy_for_implicit_arbitrary_stepper_t<
    state_t, jacobian_t, system_t, scalar_t, Args...>;
  using jacobian_policy_t = ::pressio::mpl::variadic::at_or_t<void, ic4::value, Args...>;

  // make sure we have non-void policies
  static_assert( !std::is_void<residual_policy_t>::value and
		 !std::is_void<jacobian_policy_t>::value,
  		 "If you want to use an arbitrary implicit stepper, \
i.e. when ode::ImplicitEnum::Arbitrary, you must provide user-defined \
residual and jacobian policies which are in charge of computing the \
time discrete residual and its jacobian. ");

  // now that the policies are set, find order and num of aux states needed.
  // detect order and num aux states from the residual policy
  static constexpr types::stepper_order_t order_value_from_res_pol     = residual_policy_t::stepper_order;
  static constexpr types::stepper_n_states_t n_aux_states_from_res_pol = residual_policy_t::num_aux_states;
  // detect order and num aux states from the jacobian policy
  static constexpr types::stepper_order_t order_value_from_jac_pol     = jacobian_policy_t::stepper_order;
  static constexpr types::stepper_n_states_t n_aux_states_from_jac_pol = jacobian_policy_t::num_aux_states;

  // make sure they are the same otherwise it might be that residual and jacobian are
  // doing different things
  static_assert( order_value_from_res_pol == order_value_from_jac_pol,
		 "The stepper order I detected from the residual policy does not match the \
one detected from the jacobian policy. This can be because you are using two \
policies together that actually do things differently. So one of them has to be fixed because \
the stepper order should match." );

  static_assert( n_aux_states_from_res_pol == n_aux_states_from_jac_pol,
		 "The num_aux_states I detected from the residual policy does not match the \
one detected from the jacobian policy. This can be because you are using two \
policies together that actually do things differently. So one of them has to be fixed because \
the num of aux states neeeded should match." );

  // if we get here, we can use order and num states from either policies since they are the same
  static constexpr auto order_value = order_value_from_res_pol;
  static constexpr auto numAuxStates = n_aux_states_from_res_pol;

//   // we need to make sure that the residual and jac policies are not standard policies
//   using policy_picker = impl::StdPoliciesPicker<system_t, state_t, residual_t, jacobian_t>;
//   using standard_res_policy_t = typename policy_picker::standard_res_policy_t;
//   using standard_jac_policy_t = typename policy_picker::standard_jac_policy_t;
//   static_assert( mpl::not_same<residual_policy_t, residual_policy_t>::value and
// 		 mpl::not_same<jacobian_policy_t, jacobian_policy_t>::value ,
//   		 "You cannot use standard policies to an arbitrary implicit stepper \
// i.e. when ode::ImplicitEnum::Arbitrary, you must provide user-defined \
// residual and jacobian policies which are in charge of computing the \
// time discrete residual and its jacobian. ");

};

}}}//end namespace pressio::ode::details
#endif
