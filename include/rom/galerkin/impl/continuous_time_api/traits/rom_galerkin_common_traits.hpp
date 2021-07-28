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

#ifndef ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_TRAITS_ROM_GALERKIN_COMMON_TRAITS_HPP_
#define ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_TRAITS_ROM_GALERKIN_COMMON_TRAITS_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template<typename tag>
struct _num_fom_states_needed{
  static constexpr std::size_t value = 1;
};

template<>
struct _num_fom_states_needed<::pressio::ode::implicitmethods::Euler>{
  static constexpr std::size_t value = 1;
};

template<>
struct _num_fom_states_needed<::pressio::ode::implicitmethods::BDF2>{
  static constexpr std::size_t value = 1;
};

template<>
struct _num_fom_states_needed<::pressio::ode::implicitmethods::CrankNicolson>{
  static constexpr std::size_t value = 2;
};

template <
  typename ode_tag,
  typename fom_system_type,
  typename rom_state_type,
  typename decoder_type,
  typename ud_ops_type
  >
struct CommonTraitsContinuousTimeApi
{
  using fom_system_t	= fom_system_type;
  using scalar_t	= typename fom_system_t::scalar_type;

  // rom state and native type
  using galerkin_state_t	= rom_state_type;
  using galerkin_native_state_t	=
    typename ::pressio::containers::details::traits<galerkin_state_t>::wrapped_t;

  // ---------------------
  // verify decoder
  static_assert
  (::pressio::rom::constraints::decoder<decoder_type, galerkin_state_t>::value,
   "A valid decoder type must be passed to define a Galerkin problem");
  using decoder_t = decoder_type;
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
  // now we don't allow fom state and fom velocity to have different types
  // but need to make sure this assumption is consistent with fom class
  using fom_velocity_t = fom_state_t;
  using fom_native_velocity_t = typename fom_system_type::velocity_type;
  static_assert
  (std::is_same<fom_native_state_t, fom_native_velocity_t>::value,
   "Currently, the fom velocity type must be the same as the state type.");

  // ---------------------
  // fom state reconstructor type
  using fom_state_reconstr_t =
    typename ::pressio::rom::impl::FomStateReconHelper<
    ud_ops_type>::template type<scalar_t, fom_state_t, decoder_t>;

  // ---------------------------------------------------------------
  /* the fom states manager
     For explicit time stepping, we only need to store one fom state that
     we reconstruct every time we need to compute the FOM velocity
     Keep in mind that this is Galerkin, so the FOM states are only needed
     when we query the FOM velocity.
     So for implicit time stepping, we might need to store more than one
     depending on how many FOM velocity evaluations are needed.
     For BDF1 and BDF2, we only need one evaluation of the
     FOM velocity at n+1, so we only need one FOM state.
     For CrankNicolson, we need two evaluations of the FOM velocity
     at n and n+1, so we need two FOM states.
   */
  static constexpr auto nstates = _num_fom_states_needed<ode_tag>::value;
  using tagtype = mpl::conditional_t<
    ::pressio::ode::predicates::is_explicit_stepper_tag<ode_tag>::value,
    ::pressio::rom::UnsteadyExplicit, ::pressio::rom::UnsteadyImplicit>;
  using fom_states_manager_t = ::pressio::rom::ManagerFomStates<
    tagtype, fom_state_t, fom_state_reconstr_t, ud_ops_type, nstates>;

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

}}}}//end  namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_TRAITS_ROM_GALERKIN_COMMON_TRAITS_HPP_
