/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_type_generator_common_velocity_api.hpp
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

#ifndef ROM_GALERKIN_TYPE_GENERATOR_COMMON_VELOCITY_API_HPP_
#define ROM_GALERKIN_TYPE_GENERATOR_COMMON_VELOCITY_API_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <typename fom_t, typename galerkin_state_t, typename enable = void>
struct ExtractNativeHelper;

template <typename fom_t, typename galerkin_state_t>
struct ExtractNativeHelper<
  fom_t, galerkin_state_t
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  ,mpl::enable_if_t<
    and !::pressio::containers::predicates::is_vector_wrapper_pybind<galerkin_state_t>::value
    >
#endif
  >
{
  using fom_native_state_t    = typename fom_t::state_type;
  using fom_native_velocity_t = typename fom_t::velocity_type;
};

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <typename fom_t, typename galerkin_state_t>
struct ExtractNativeHelper<
  fom_t, galerkin_state_t,
  mpl::enable_if_t<
    ::pressio::containers::predicates::is_vector_wrapper_pybind<galerkin_state_t>::value
    >
  >
{
  using fom_native_state_t    = typename ::pressio::containers::details::traits<galerkin_state_t>::wrapped_t;
  using fom_native_velocity_t = typename ::pressio::containers::details::traits<galerkin_state_t>::wrapped_t;
};
#endif
// ------------------------------------------------------------------------------------

template <typename ops_t, typename enable = void>
struct FomStateReconHelper;

template <typename ops_t>
struct FomStateReconHelper<
  ops_t, mpl::enable_if_t< std::is_void<ops_t>::value >
  >
{
  template <typename scalar_t, typename fom_state_t, typename decoder_t>
  using type = FomStateReconstructor<scalar_t, fom_state_t, decoder_t>;
};

template <typename ops_t>
struct FomStateReconHelper<
  ops_t, mpl::enable_if_t< !std::is_void<ops_t>::value >
  >
{
  template <typename scalar_t, typename fom_state_t, typename decoder_t>
  using type = FomStateReconstructor<scalar_t, fom_state_t, decoder_t, ops_t>;
};
//------------------------------------------------------------------------------


template <typename fom_type, typename rom_state_type, typename ...Args >
struct CommonTraitsContinuousTimeApi
{
  // the scalar type
  using scalar_t		= typename ::pressio::containers::details::traits<rom_state_type>::scalar_t;

  using fom_t			= fom_type;
  using fom_native_state_t	= typename ExtractNativeHelper<fom_t, rom_state_type>::fom_native_state_t;
  using fom_native_velocity_t	= typename ExtractNativeHelper<fom_t, rom_state_type>::fom_native_velocity_t;
  // fom wrapper types
  using fom_state_t		= ::pressio::containers::Vector<fom_native_state_t>;
  using fom_velocity_t		= ::pressio::containers::Vector<fom_native_velocity_t>;

  // rom state type and native type
  using galerkin_state_t	= rom_state_type;
  using galerkin_native_state_t	= typename ::pressio::containers::details::traits<galerkin_state_t>::wrapped_t;

  // the Galerkin rhs type is (for now) same as state type
  using galerkin_residual_t	= galerkin_state_t;

  // verify the sequence contains a valid decoder type
  using ic2 = ::pressio::mpl::variadic::find_if_ternary_pred_t<
    galerkin_state_t, fom_state_t, ::pressio::rom::concepts::admissible_decoder, Args...>;
  using decoder_t = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert(!std::is_void<decoder_t>::value and ic2::value < sizeof... (Args),
		"A valid decoder type must be passed to define a ROM Galerkin problem");
  using decoder_jac_t = typename decoder_t::jacobian_type;

  // if we have an admissible user-defined ops
  using icUdOps = ::pressio::mpl::variadic::find_if_quaternary_pred_t<
    decoder_jac_t, galerkin_state_t, fom_state_t,
    ::pressio::rom::concepts::custom_ops_galerkin_continuous_time, Args...>;
  using ud_ops_t = ::pressio::mpl::variadic::at_or_t<void, icUdOps::value, Args...>;

  // fom state reconstructor type
  using fom_state_reconstr_t =
    typename FomStateReconHelper<ud_ops_t>::template type<scalar_t, fom_state_t, decoder_t>;

  // class type holding fom states data
  using fom_states_manager_t = ::pressio::rom::ManagerFomStatesStatic<fom_state_t, 1, fom_state_reconstr_t, ud_ops_t>;
};

}}}}//end  namespace pressio::rom::galerkin::impl
#endif
