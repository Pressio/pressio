/*
//@HEADER
// ************************************************************************
//
// rom_compose_lspg_impl.hpp
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

#ifndef ROM_LSPG_IMPL_ROM_COMPOSE_STEADY_LSPG_IMPL_HPP_
#define ROM_LSPG_IMPL_ROM_COMPOSE_STEADY_LSPG_IMPL_HPP_

#include "../rom_lspg_problem_tags.hpp"

#include "./policies/rom_lspg_steady_residual_policy.hpp"
#include "./policies/rom_lspg_steady_jacobian_policy.hpp"

#include "rom_lspg_steady_system.hpp"

#include "./traits/rom_lspg_steady_common_traits.hpp"
#include "./traits/rom_lspg_steady_default_problem_traits.hpp"
#include "./traits/rom_lspg_steady_preconditioned_problem_traits.hpp"

#include "./rom_lspg_steady_default_problem.hpp"
#include "./rom_lspg_steady_preconditioned_problem.hpp"


namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template<typename problem_tag, typename enable, typename fom_system_type, typename ...Args>
struct composeSteady
{
  //if we are here, something is wrong, find out what
  static_assert
  (::pressio::rom::find_discrepancies_with_steady_system_api<fom_system_type>::value, "");
  using type = void;
};

//------------------------
// specialize compose
//------------------------

// default
template<
  typename fom_system_type,
  typename lspg_state_type,
  typename decoder_type
  >
struct composeSteady<
  ::pressio::rom::lspg::impl::Default,
  mpl::enable_if_t<
    ::pressio::rom::concepts::steady_system<fom_system_type>::value
    >,
  fom_system_type, lspg_state_type, decoder_type>
{
  using type = ::pressio::rom::lspg::impl::steady::DefaultProblemSteady<
      fom_system_type, lspg_state_type, decoder_type>;
};

// precond
template<
  typename fom_system_type,
  typename lspg_state_type,
  typename decoder_type,
  typename precond_type
  >
struct composeSteady<
  ::pressio::rom::lspg::impl::Preconditioned,
  mpl::enable_if_t<
    ::pressio::rom::concepts::steady_system<fom_system_type>::value
    >,
  fom_system_type, lspg_state_type, decoder_type, precond_type>
{
  using type = ::pressio::rom::lspg::impl::steady::PreconditionedProblemSteady<
      fom_system_type, lspg_state_type, decoder_type, precond_type>;
};

}}}}
#endif  // ROM_LSPG_IMPL_ROM_COMPOSE_LSPG_IMPL_HPP_
