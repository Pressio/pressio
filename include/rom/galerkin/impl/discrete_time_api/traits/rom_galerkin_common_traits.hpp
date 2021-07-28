/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_common_traits.hpp
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

#ifndef ROM_GALERKIN_IMPL_DISCRETE_TIME_API_TRAITS_ROM_GALERKIN_COMMON_TRAITS_HPP_
#define ROM_GALERKIN_IMPL_DISCRETE_TIME_API_TRAITS_ROM_GALERKIN_COMMON_TRAITS_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <
  typename fom_system_type,
  typename rom_state_type,
  typename decoder_type,
  typename ...Args
  >
struct CommonTraitsDiscreteTimeApi
{
  using fom_system_t		= fom_system_type;
  using scalar_t		= typename fom_system_t::scalar_type;

  // rom types
  using galerkin_state_t	= rom_state_type;
  using galerkin_native_state_t	=
    typename ::pressio::containers::details::traits<galerkin_state_t>::wrapped_t;

  // ---------------------
  // verify decoder
  static_assert
  (::pressio::rom::constraints::decoder<decoder_type, galerkin_state_t>::value,
   "A valid decoder type must be passed to define a Galerkin problem");
  using decoder_t     = decoder_type;
  using decoder_jac_t = typename decoder_type::jacobian_type;

  // ---------------------
  // detect fom state type (supposed to be wrapper) from decoder
  // ensure it is consistent with the (native) fom_state_type from the app
  using fom_state_t = typename decoder_type::fom_state_type;
  using fom_native_state_t = typename fom_system_type::state_type;
  static_assert
  (std::is_same<
   typename ::pressio::containers::details::traits<fom_state_t>::wrapped_t,
   fom_native_state_t>::value,
   "The fom state type detected in the fom class must match the fom state type used in the decoder");

  // ---------------------
  // for now we don't allow state and residual to have different types
  // but need to make sure this assumption is consistent with fom class
  using fom_residual_t		= fom_state_t;
  using fom_native_residual_t = typename fom_system_type::discrete_time_residual_type;
  static_assert
  (std::is_same<fom_native_state_t, fom_native_residual_t>::value,
   "Currently, the fom discrete time residual type must be the same as the state type.");

  // ---------------------
  /* fom_apply_jacobian_t is type of J*decoder_jac_t where
   * * J is the jacobian of the fom rhs
   * * decoder_jac_t is the type of the decoder jacobian
   * For now, set fom_apply_jacobian_t to be of same type as decoder_jac_t
   * not bad assumption since all matrices are left-applied to decoder_jac_t
   */
  using fom_apply_jacobian_t = decoder_jac_t;

  // if we have an admissible user-defined ops
  using ud_ops_t = void;

  // ---------------------
  // fom state reconstructor type
  using fom_state_reconstr_t = FomStateReconstructor<scalar_t, fom_state_t, decoder_t>;

  //-------------------------------
  // find the order setter in Args
  //-------------------------------
  using ic1 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::ode::predicates::IsStepperOrderSetter, Args...>;
  using order_setter = ::pressio::mpl::variadic::at_or_t<void, ic1::value, Args...>;
  static_assert( !std::is_void<order_setter>::value,
  		 "To use Galerkin with residual api, you need to set the order of the stepper \n \
at compile time by passing to a template argument as follows: \n \
::pressio::ode::types::StepperOrder<your_order_value>.");
  // store
  static constexpr ::pressio::ode::types::stepper_order_t order_value = order_setter::value;

  //-----------------------------------------------------------
  // find the total number of states needed
  //-----------------------------------------------------------
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::ode::predicates::IsStepperTotalNumStatesSetter, Args...>;
  using tot_n_setter = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert( !std::is_void<tot_n_setter>::value,
  		 "\nTo use Galerkin with residual api, you need to set the \
total number of states needed for the stepper at compile time by passing \
to a template argument as follows: \n \
::pressio::ode::types::StepperTotalNumberOfStates<your_order_value>. \n \
Note that this is the total number of states needed including previous ones, \n \
basically the size of the stpper stencil.");

  // total number of fom states needed (size of stencil plus the state at current step)
  static constexpr std::size_t numstates = tot_n_setter::value;
  using fom_states_manager_t = ::pressio::rom::ManagerFomStates<
    ::pressio::rom::UnsteadyImplicit, fom_state_t, fom_state_reconstr_t,
    ud_ops_t, numstates>;

  // ---------------------
  // sentinel to tell if we are doing bindings for p4py:
  // always false if pybind is disabled, otherwise detect from galerkin state
  static constexpr bool binding_sentinel =
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    ::pressio::containers::predicates::is_tensor_wrapper_pybind<galerkin_state_t>::value;
#else
  false;
#endif
};

}}}}//end  namespace
#endif  // ROM_GALERKIN_IMPL_DISCRETE_TIME_API_TRAITS_ROM_GALERKIN_COMMON_TRAITS_HPP_
