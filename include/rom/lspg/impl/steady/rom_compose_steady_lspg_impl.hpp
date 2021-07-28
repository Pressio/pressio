/*
//@HEADER
// ************************************************************************
//
// rom_compose_steady_lspg_impl.hpp
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

#ifndef ROM_LSPG_IMPL_STEADY_ROM_COMPOSE_STEADY_LSPG_IMPL_HPP_
#define ROM_LSPG_IMPL_STEADY_ROM_COMPOSE_STEADY_LSPG_IMPL_HPP_

#include "../rom_problem_tags.hpp"

#include "./policies/rom_lspg_steady_residual_policy.hpp"
#include "./policies/rom_lspg_steady_jacobian_policy.hpp"

#include "rom_lspg_steady_system.hpp"

#include "./traits/rom_lspg_steady_common_traits.hpp"
#include "./traits/rom_lspg_steady_default_problem_traits.hpp"
#include "./traits/rom_lspg_steady_preconditioned_problem_traits.hpp"
#include "./traits/rom_lspg_steady_masked_problem_traits.hpp"
#include "./traits/rom_lspg_steady_hyper_reduced_problem_traits.hpp"
#include "./traits/rom_lspg_steady_preconditioned_hyper_reduced_problem_traits.hpp"

#include "./rom_lspg_steady_default_problem.hpp"
#include "./rom_lspg_steady_preconditioned_problem.hpp"
#include "./rom_lspg_steady_masked_problem.hpp"
#include "./rom_lspg_steady_hyper_reduced_problem.hpp"
#include "./rom_lspg_steady_preconditioned_hyper_reduced_problem.hpp"

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template<
  class problem_tag,
  class enable,
  class fom_system_type,
  class decoder_type,
  class ...Args
  >
struct composeSteady
{
  //if we are here, something is wrong, find out what
  static_assert
  (::pressio::rom::why_not_steady_system_with_user_provided_apply_jacobian
   <fom_system_type, typename decoder_type::jacobian_type
   >::value, "");

  using type = void;
};

//------------------------
// specialize compose
//------------------------

// default
template<
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type
  >
struct composeSteady<
  ::pressio::rom::lspg::impl::Default,
  mpl::enable_if_t<
    ::pressio::rom::constraints::steady_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  fom_system_type, decoder_type, lspg_state_type>
{
  using type = ::pressio::rom::lspg::impl::steady::DefaultProblemSteady<
    fom_system_type, lspg_state_type, decoder_type>;
};

// preconditioned default
template<
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type,
  typename precond_type
  >
struct composeSteady<
  ::pressio::rom::lspg::impl::Preconditioned,
  mpl::enable_if_t<
    ::pressio::rom::constraints::steady_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  fom_system_type, decoder_type, lspg_state_type, precond_type>
{
  using type = ::pressio::rom::lspg::impl::steady::PreconditionedProblemSteady<
    fom_system_type, lspg_state_type, decoder_type, precond_type>;
};

// masked
template<
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type,
  typename masker_type
  >
struct composeSteady<
  ::pressio::rom::lspg::impl::Masked,
  mpl::enable_if_t<
    ::pressio::rom::constraints::steady_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  fom_system_type, decoder_type, lspg_state_type, masker_type>
{
  using type = ::pressio::rom::lspg::impl::steady::MaskedProblemSteady<
    fom_system_type, lspg_state_type, decoder_type, masker_type>;
};

// hyper-reduced
template<
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type
  >
struct composeSteady<
  ::pressio::rom::lspg::impl::HyperReduced,
  mpl::enable_if_t<
    ::pressio::rom::constraints::steady_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  fom_system_type, decoder_type, lspg_state_type
  >
{
  using type = ::pressio::rom::lspg::impl::steady::HyperReducedProblemSteady<
    fom_system_type, lspg_state_type, decoder_type>;
};

// preconditioned hyper-reduced
template<
  typename fom_system_type,
  typename decoder_type,
  typename lspg_state_type,
  typename precond_type
  >
struct composeSteady<
  ::pressio::rom::lspg::impl::PreconditionedHyperReduced,
  mpl::enable_if_t<
    ::pressio::rom::constraints::steady_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  fom_system_type, decoder_type, lspg_state_type, precond_type>
{
  using type = ::pressio::rom::lspg::impl::steady::PreconditionedHyperReducedProblemSteady<
    fom_system_type, lspg_state_type, decoder_type, precond_type>;
};

}}}}
#endif  // ROM_LSPG_IMPL_STEADY_ROM_COMPOSE_STEADY_LSPG_IMPL_HPP_
