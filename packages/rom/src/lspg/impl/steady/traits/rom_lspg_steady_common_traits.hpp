/*
//@HEADER
// ************************************************************************
//
// rom_lspg_steady_common_traits.hpp
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

#ifndef ROM_LSPG_IMPL_STEADY_TRAITS_ROM_LSPG_STEADY_COMMON_TRAITS_HPP_
#define ROM_LSPG_IMPL_STEADY_TRAITS_ROM_LSPG_STEADY_COMMON_TRAITS_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{ namespace steady{

template <
  typename fom_system_type,
  typename lspg_state_type,
  typename decoder_type
  >
struct CommonTraits
{
  static_assert
  (::pressio::rom::lspg::constraints::state<lspg_state_type>::value,
   "The lspg_state_type is not a valid rom state");

  using fom_system_t	      = fom_system_type;
  using scalar_t              = typename fom_system_t::scalar_type;

  // lspg state types
  using lspg_state_t	    = lspg_state_type;
  using lspg_native_state_t = typename ::pressio::containers::details::traits<lspg_state_t>::wrapped_t;

  // ---------------------
  // check for valid decoder
  static_assert
  (::pressio::rom::constraints::decoder<decoder_type, lspg_state_t>::value,
   "A valid decoder type must be passed to define a LSPG problem");
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
  using fom_native_residual_t = typename fom_system_type::residual_type;
  static_assert
  (std::is_same<
   typename ::pressio::containers::details::traits<fom_residual_t>::wrapped_t,
   fom_native_residual_t>::value,
   "Currently, the fom discrete time residual type must be the same as the state type.");

  // ---------------------
  /* for steady lspg, the lspg residual is literally the rhs (aka velocity fo us) */
  using lspg_residual_t = fom_residual_t;

  /* lspg_jacobian_t is type to represent dR/dx_rom where R is the residual
   * dR/dx_rom = df/dxFom * dxFom/dxRom
   *  so df/dxFom = J  and dxFom/dxRom = decoder_jacobian
   * For now, set lspg_jacobian_t to be of same type as decoder_jac_t
   * not a bad assumption since all matrices are left-applied to decoder_jac_t
   */
  using lspg_jacobian_t	= decoder_jac_t;

  // ---------------------
  // fro now, later on this needs to be detected from args
  using ud_ops_t = void;

  // ---------------------
  // fom state reconstructor type
  using fom_state_reconstr_t =
    FomStateReconstructor<scalar_t, fom_state_t, decoder_t>;

  // for steady lspg, we only need to store one FOM state
  using fom_states_manager_t = ::pressio::rom::ManagerFomStates<
    ::pressio::rom::Steady, fom_state_t, fom_state_reconstr_t, ud_ops_t>;

  // sentinel to tell if we are doing bindings for p4py:
  // always false if pybind is disabled, otherwise detect from rom state
  static constexpr bool binding_sentinel =
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    ::pressio::containers::predicates::is_tensor_wrapper_pybind<lspg_state_t>::value;
#else
  false;
#endif
};

}}}}}
#endif  // ROM_LSPG_IMPL_STEADY_TRAITS_ROM_LSPG_STEADY_COMMON_TRAITS_HPP_
