/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_common_traits_velocity_api.hpp
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

#ifndef ROM_LSPG_UNSTEADY_COMMON_TRAITS_VELOCITY_API_HPP_
#define ROM_LSPG_UNSTEADY_COMMON_TRAITS_VELOCITY_API_HPP_

#include "../../shared/rom_lspg_unsteady_aux_stepper_type_helper.hpp"
#include "../../shared/rom_lspg_unsteady_fom_states_storage_capacity_helper.hpp"
#include "../../shared/rom_lspg_unsteady_fom_state_reconstructor_helper.hpp"

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace unsteady{

template <typename system_t, typename lspg_state_t, typename enable = void>
struct ExtractNativeHelper;

template <typename system_t, typename lspg_state_t>
struct ExtractNativeHelper<
  system_t, lspg_state_t
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  ,mpl::enable_if_t<
     !::pressio::containers::predicates::is_vector_wrapper_pybind<lspg_state_t>::value
    >
#endif
  >
{
  using fom_native_state_t    = typename system_t::state_type;
  using fom_native_velocity_t = typename system_t::velocity_type;
};

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <typename system_t, typename lspg_state_t>
struct ExtractNativeHelper<
  system_t, lspg_state_t,
  mpl::enable_if_t<
    ::pressio::containers::predicates::is_vector_wrapper_pybind<lspg_state_t>::value
    >
  >
{
  using fom_native_state_t    = typename ::pressio::containers::details::traits<lspg_state_t>::wrapped_t;
  using fom_native_velocity_t = typename ::pressio::containers::details::traits<lspg_state_t>::wrapped_t;
};
#endif
// ------------------------------------------------------------------------------------


template <
  typename stepper_tag,
  typename system_type,
  typename lspg_state_type,
  typename ...Args>
struct CommonTraitsContinuousTimeApi
{
  // the scalar type
  using scalar_t = typename ::pressio::containers::details::traits<lspg_state_type>::scalar_t;

  using system_t			= system_type;
  using fom_native_state_t	= typename ExtractNativeHelper<system_t, lspg_state_type>::fom_native_state_t;
  using fom_native_velocity_t	= typename ExtractNativeHelper<system_t, lspg_state_type>::fom_native_velocity_t;

  // fom wrapper types
  using fom_state_t		= ::pressio::containers::Vector<fom_native_state_t>;
  using fom_velocity_t		= ::pressio::containers::Vector<fom_native_velocity_t>;

  // rom state type (passed in)
  using lspg_state_t		= lspg_state_type;

  /* for LSPG, the lspg_residual_type = fom_velocity_type */
  using lspg_residual_t		= fom_velocity_t;

  // verify that args contains a valid decoder type
  using ic2 = ::pressio::mpl::variadic::find_if_ternary_pred_t<
    lspg_state_t, fom_state_t, ::pressio::rom::concepts::admissible_decoder, Args...>;
  using decoder_t = ::pressio::mpl::variadic::at_or_t<void, ic2::value, Args...>;
  static_assert(!std::is_void<decoder_t>::value and ic2::value < sizeof... (Args),
		"A valid decoder type must be passed to define a LSPG problem");
  using decoder_jac_t = typename decoder_t::jacobian_type;

  /* lspg_matrix_t is type of J*decoder_jac_t (in the most basic case) where
   * * J is the jacobian of the fom rhs
   * * decoder_jac_t is the type of the decoder jacobian
   * In more complex cases, we might have (something)*J*decoder_jac_t,
   * where (something) is product of few matrices.
   * For now, set lspg_matrix_t to be of same type as decoder_jac_t
   * if phi is MV<>, then lspg_matrix_t = containers::MV<>
   * if phi is Matrix<>, then we have containers::Matrix<>
   * not bad assumption since all matrices are left-applied to decoder_jac_t
   */
  using lspg_matrix_t		= decoder_jac_t;

  // if we have an admissible user-defined ops
  using icUdOps = ::pressio::mpl::variadic::find_if_quaternary_pred_t<
    decoder_jac_t, lspg_state_t, fom_state_t,
    ::pressio::rom::concepts::custom_ops_for_unsteady_lspg_velocity_api, Args...>;
  using ud_ops_t = ::pressio::mpl::variadic::at_or_t<void, icUdOps::value, Args...>;

  // fom state reconstructor type
  using fom_state_reconstr_t =
    typename FomStateReconHelper<ud_ops_t>::template type<scalar_t, fom_state_t, decoder_t>;

  // total num of fom states (i.e. stencil size plus the state at current step) needed
  static constexpr auto numStates = fomStatesStorageCapacityHelper<stepper_tag>::value;

  // type of class holding the fom states
  using fom_states_manager_t = ::pressio::rom::ManagerFomStatesStatic<fom_state_t, numStates, fom_state_reconstr_t, ud_ops_t>;
};

}}}}}//end  namespace pressio::rom::lspg::unstedy::impl
#endif
