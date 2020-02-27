/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_type_generator_common_impl.hpp
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

#ifndef ROM_GALERKIN_TYPE_GENERATOR_COMMON_IMPL_HPP_
#define ROM_GALERKIN_TYPE_GENERATOR_COMMON_IMPL_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <typename fom_t, typename galerkin_state_t, typename enable = void>
struct ExtractNativeHelper;

template <typename fom_t, typename galerkin_state_t>
struct ExtractNativeHelper<
  fom_t, galerkin_state_t,
  mpl::enable_if_t<
    ::pressio::rom::meta::is_legitimate_model_for_galerkin<fom_t>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    and !::pressio::containers::meta::is_vector_wrapper_pybind<galerkin_state_t>::value
#endif
    >
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
    ::pressio::rom::meta::is_legitimate_model_for_galerkin<fom_t>::value
    and ::pressio::containers::meta::is_vector_wrapper_pybind<galerkin_state_t>::value
    >
  >
{
  using fom_native_state_t    = typename ::pressio::containers::details::traits<galerkin_state_t>::wrapped_t;
  using fom_native_velocity_t = typename ::pressio::containers::details::traits<galerkin_state_t>::wrapped_t;
};
#endif
// ------------------------------------------------------------------------------------

template < typename galerkin_state_type, typename ...Args >
struct GalerkinCommonTypes
{
  // rom state type
  using galerkin_state_t	= galerkin_state_type;
  // native type
  using galerkin_native_state_t	= typename ::pressio::containers::details::traits<galerkin_state_t>::wrapped_t;

  // verify the sequence contains a valid decoder type
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::rom::meta::is_legitimate_decoder_type, Args...>;
  using decoder_t = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert(!std::is_void<decoder_t>::value and ic2::value < sizeof... (Args),
		"A valid decoder type must be passed to define a ROM Galerkin problem");
  using decoder_jac_t = typename decoder_t::jacobian_t;

  // verify that args contains a valid fom/adapter type for Galerkin
  using ic1 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::rom::meta::is_legitimate_model_for_galerkin, Args...>;
  using fom_t = ::pressio::mpl::variadic::at_or_t<void, ic1::value, Args...>;
  static_assert(!std::is_void<fom_t>::value and ic1::value < sizeof... (Args),
		"A valid adapter/fom type must be passed to define a ROM Galerkin problem");

  // the scalar type
  using scalar_t = typename ::pressio::containers::details::traits<galerkin_state_type>::scalar_t;

  // the GALERKIN rhs type is (for now) same as state type
  using galerkin_residual_t	= galerkin_state_t;

  // get the native types
  using fom_native_state_t	= typename ExtractNativeHelper<fom_t, galerkin_state_t>::fom_native_state_t;
  using fom_native_velocity_t	= typename ExtractNativeHelper<fom_t, galerkin_state_t>::fom_native_velocity_t;
  // declare wrapper types
  using fom_state_t		= ::pressio::containers::Vector<fom_native_state_t>;
  using fom_velocity_t		= ::pressio::containers::Vector<fom_native_velocity_t>;

  // fom state reconstructor type
  using fom_state_reconstr_t	= FomStateReconstructor<scalar_t, fom_state_t, decoder_t>;

  // class type holding fom states data
  using fom_states_data = ::pressio::rom::FomStatesStaticContainer<fom_state_t, 1, fom_state_reconstr_t>;

  // for now, set ops to void, i.e. we only use pressio ops
  using ud_ops_t = void;
};

}}}}//end  namespace pressio::rom::galerkin::impl
#endif
