/*
//@HEADER
// ************************************************************************
//
// rom_compose_lspg.hpp
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

#ifndef ROM_LSPG_ROM_COMPOSE_LSPG_HPP_
#define ROM_LSPG_ROM_COMPOSE_LSPG_HPP_

#include "./impl/steady/rom_compose_steady_lspg_impl.hpp"
#include "./impl/unsteady/rom_compose_unsteady_lspg_impl.hpp"

namespace pressio{ namespace rom{ namespace lspg{

/*default

  unsteady:
  template<stepper_tag, fom_type, romstate_t, decoder_t>
  template<stepper_tag, fom_type, romstate_t, decoder_t, udops_t>

  steady:
  template<fom_type, romstate_t, decoder_t>
*/
template<typename T1, typename ...Args>
using composeDefaultProblem =
  typename std::conditional<
  ::pressio::ode::predicates::is_stepper_tag<T1>::value,
  impl::composeUnsteady<impl::Default, void, T1, Args...>,
  impl::composeSteady<impl::Default, void, T1, Args...>
  >::type;


/*preconditioned

  unsteady:
  template<stepper_tag, fom_type, romstate_t, decoder_t, precond_t>

  cases for steady:
  template<fom_type, romstate_t, decoder_t, precond_t>
*/
template<typename T1, typename ...Args>
using composePreconditionedProblem =
  typename std::conditional<
  ::pressio::ode::predicates::is_stepper_tag<T1>::value,
  impl::composeUnsteady<impl::Preconditioned, void, T1, Args...>,
  impl::composeSteady<impl::Preconditioned, void, T1, Args...>
  >::type;


/*masked

  unsteady:
  template<stepper_tag, fom_type, romstate_t, decoder_t, masker_t>

  steady: TBI
*/
template<typename T1, typename ...Args>
using composeMaskedProblem =
  typename std::conditional<
  ::pressio::ode::predicates::is_stepper_tag<T1>::value,
  impl::composeUnsteady<impl::Masked, void, T1, Args...>,
  void
  >::type;


// /*hyperreduced
//   for now just an alias to above
// */
// template<typename T1, typename ...Args>
// using composeHyperreducedProblem = composeDefaultProblem<T1, Args...>;

// /*preconditioned hyperreduced */
// template<typename T1, typename ...Args>
// using composePreconditionedHyperreducedProblem =
//   composePreconditionedProblem<T1, Args...>;

}}}
#endif  // ROM_LSPG_ROM_COMPOSE_LSPG_HPP_
