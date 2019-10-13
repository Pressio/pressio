/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_type_generator_common.hpp
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

#ifndef ROM_GALERKIN_TYPE_GENERATOR_COMMON_HPP_
#define ROM_GALERKIN_TYPE_GENERATOR_COMMON_HPP_

#include "../../rom_ConfigDefs.hpp"
#include "../../rom_fwd.hpp"
#include "../../rom_data_fom_rhs.hpp"
#include "../../rom_data_fom_states.hpp"
#include "../../policies/rom_evaluate_fom_velocity_unsteady_policy.hpp"
#include "../../policies/rom_apply_fom_jacobian_unsteady_policy.hpp"
#include "../../../../ode/src/ode_fwd.hpp"
#include "../../meta/rom_is_legitimate_model_for_galerkin.hpp"
#include "../../meta/rom_is_legitimate_decoder_type.hpp"

namespace pressio{ namespace rom{ namespace impl{

template < bool doingPython, typename galerkin_state_type, typename ...Args >
struct GalerkinCommonTypes;

template < typename galerkin_state_type, typename ...Args >
struct GalerkinCommonTypes<false, galerkin_state_type, Args...>
{
  // verify that args contains a valid fom/adapter type for Galerkin
  using ic1 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::rom::meta::is_legitimate_model_for_galerkin, Args...>;
  using fom_t = ::pressio::mpl::variadic::at_or_t<void, ic1::value, Args...>;
  static_assert(!std::is_void<fom_t>::value and ic1::value < sizeof... (Args),
		"A valid adapter/fom type must be passed to define a ROM Galerkin problem");

  // get the native types from the full-order model (fom)
  using scalar_t		= typename fom_t::scalar_type;
  using fom_native_state_t	= typename fom_t::state_type;
  using fom_native_velocity_t	= typename fom_t::velocity_type;

  // declare fom wrapper types
  using fom_state_t		= ::pressio::containers::Vector<fom_native_state_t>;
  using fom_velocity_t		= ::pressio::containers::Vector<fom_native_velocity_t>;

  // rom state type (passed in)
  using galerkin_state_t	= galerkin_state_type;

  // the GALERKIN rhs type is (for now) same as state type
  using galerkin_residual_t	= galerkin_state_type;

  // verify the sequence contains a valid decoder type
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::rom::meta::is_legitimate_decoder_type, Args...>;
  using decoder_t = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert(!std::is_void<decoder_t>::value and ic2::value < sizeof... (Args),
		"A valid decoder type must be passed to define a ROM Galerkin problem");
  using decoder_jac_t		= typename decoder_t::jacobian_t;

  // fom state reconstructor type
  using fom_state_reconstr_t	= FomStateReconstructor<fom_state_t, decoder_t>;

  // class type holding fom states data
  using fom_states_data = ::pressio::rom::FomStatesData<
	fom_state_t, 0, fom_state_reconstr_t>;

  // class type holding fom rhs data
  using fom_velocity_data = ::pressio::rom::FomRhsData<fom_velocity_t>;

  // for now, set ops to void, i.e. we only use pressio ops
  using ud_ops_t = void;
};


#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template < typename galerkin_state_type, typename ...Args >
struct GalerkinCommonTypes<true, galerkin_state_type, Args...>
{
  // verify that args contains a valid fom/adapter type for Galerkin
  using ic1 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::rom::meta::is_legitimate_model_for_galerkin, Args...>;
  using fom_t = ::pressio::mpl::variadic::at_or_t<void, ic1::value, Args...>;
  static_assert(mpl::is_same<fom_t, pybind11::object>::value,
		"The adapter/fom type must be a pybind11::object to be valid for interfacing to Python");

  // rom state type
  using galerkin_state_t	= galerkin_state_type;

  // declare fom wrapper types
  using fom_state_t		= galerkin_state_t;
  using fom_velocity_t		= galerkin_state_t;

  // the GALERKIN residual type is (for now) same as state type
  using galerkin_residual_t	= galerkin_state_t;

  // verify the sequence contains a valid decoder type
  using ic2 = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::rom::meta::is_legitimate_decoder_type, Args...>;
  using decoder_t = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert(!std::is_void<decoder_t>::value and ic2::value < sizeof... (Args),
		"A valid decoder type must be passed to define a ROM Galerkin problem");
  using decoder_jac_t		= typename decoder_t::jacobian_t;

  // the native types
  using scalar_t		= typename decoder_t::scalar_t;
  using fom_native_state_t	= galerkin_state_t;
  using fom_native_velocity_t	= galerkin_state_t;

  // fom state reconstructor type
  using fom_state_reconstr_t	= FomStateReconstructor<fom_state_t, decoder_t>;

  // class type holding fom states data
  using fom_states_data = ::pressio::rom::FomStatesData<
	fom_state_t, 0, fom_state_reconstr_t>;

  // class type holding fom rhs data
  using fom_velocity_data = ::pressio::rom::FomRhsData<fom_velocity_t>;

  // when interfacing with Python, ops are defined by a pybind11::object
  using ud_ops_t = pybind11::object;
};
#endif

}//end namespace pressio::rom::impl


template <typename galerkin_state_type, typename ...Args>
using GalerkinCommonTypes = impl::GalerkinCommonTypes<
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  ::pressio::containers::meta::is_array_pybind11<galerkin_state_type>::value,
#else
  false,
#endif
  galerkin_state_type, Args...>;

}}//end  namespace pressio::rom
#endif
