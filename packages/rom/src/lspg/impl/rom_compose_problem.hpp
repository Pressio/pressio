/*
//@HEADER
// ************************************************************************
//
// rom_compose_problem.hpp
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

#ifndef ROM_LSPG_IMPL_ROM_COMPOSE_PROBLEM_HPP_
#define ROM_LSPG_IMPL_ROM_COMPOSE_PROBLEM_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

/*
  =====================
  === DEFAULT ===
  =====================

  steady cases:
  template<fom_type, romstate_t, decoder_t>

  unsteady cont-time api:
  template<stepper_tag, fom_type, decoder_t, romstate_t>
  template<stepper_tag, fom_type, decoder_t, romstate_t, udops_t>

  unsteady discrete-time api:
  template<stepper_tag, fom_type, decoder_t, romstate_t>
  template<stepper_tag, fom_type, decoder_t, romstate_t, udops_t>
*/

template<typename T1, typename ...Args>
using composeDefaultProblem =
  typename std::conditional<
  ::pressio::ode::predicates::is_stepper_tag<T1>::value,
  ::pressio::rom::lspg::impl::composeUnsteady<
    ::pressio::rom::lspg::impl::Default, void, T1,
    typename std::remove_cv<typename std::remove_reference<Args>::type>::type...
    >,
  ::pressio::rom::lspg::impl::composeSteady<
    ::pressio::rom::lspg::impl::Default, void, T1,
    typename std::remove_cv<typename std::remove_reference<Args>::type>::type...
    >
  >::type;

template<typename T1, typename ...Args>
using composeDefaultProblem_t = typename composeDefaultProblem<T1, Args...>::type;


/*
  ==============================
  === preconditioned-default ===
  ==============================

  for steady we support:
  template<fom_type, decoder_t, romstate_t, precond_t>

  unsteady cont-time api:
  template<stepper_tag, fom_type, decoder_t, romstate_t, precond_t>

  unsteady discrete-time api:
  template<stepper_tag, fom_type, decoder_t, romstate_t, precond_t>
*/
template<typename T1, typename ...Args>
using composePreconditionedDefaultProblem =
  typename std::conditional<
  ::pressio::ode::predicates::is_stepper_tag<T1>::value,
  ::pressio::rom::lspg::impl::composeUnsteady<
    ::pressio::rom::lspg::impl::Preconditioned, void, T1,
    typename std::remove_cv<typename std::remove_reference<Args>::type>::type...
    >,
  ::pressio::rom::lspg::impl::composeSteady<
    ::pressio::rom::lspg::impl::Preconditioned, void, T1,
    typename std::remove_cv<typename std::remove_reference<Args>::type>::type...
    >
  >::type;

template<typename T1, typename ...Args>
using composePreconditionedDefaultProblem_t =
  typename composePreconditionedDefaultProblem<T1, Args...>::type;


/*
  =====================
  === Hyper-reduced ===
  =====================

  steady:
  template<fom_type, decoder_t, romstate_t>

  unsteady continuous-time api:
  // for trilinos types we can figure out the mapping automatically
  template<stepper_tag, fom_type, decoder_t, romstate_t>
  // for shared-mem data structures where the indices are passed
  // and pressio uses these to implement hyper-reduction
  template<stepper_tag, fom_type, decoder_t, romstate_t, sample_to_stencil_t>


  unsteady discrete-time api:
  template<stepper_tag, fom_type, decoder_t, romstate_t>
*/
template<typename T1, typename ...Args>
using composeHyperReducedProblem =
  typename std::conditional<
  ::pressio::ode::predicates::is_stepper_tag<T1>::value,
  ::pressio::rom::lspg::impl::composeUnsteady<
    ::pressio::rom::lspg::impl::HyperReduced, void, T1,
    typename std::remove_cv<typename std::remove_reference<Args>::type>::type...
    >,
  ::pressio::rom::lspg::impl::composeSteady<
    ::pressio::rom::lspg::impl::HyperReduced, void, T1,
    typename std::remove_cv<typename std::remove_reference<Args>::type>::type...
    >
  >::type;

template<typename T1, typename ...Args>
using composeHyperReducedProblem_t =
  typename composeHyperReducedProblem<T1, Args...>::type;


/*
  ====================================
  === preconditioned hyper reduced ===
  ====================================

  for steady we support:
  template<fom_type, decoder_t, romstate_t, precond_t>

  unsteady cont-time api:
  template<stepper_tag, fom_type, decoder_t, romstate_t, precond_t>

  unsteady discrete-time api:
  template<stepper_tag, fom_type, decoder_t, romstate_t, precond_t>
*/
template<typename T1, typename ...Args>
using composePreconditionedHyperReducedProblem =
  typename std::conditional<
  ::pressio::ode::predicates::is_stepper_tag<T1>::value,
  ::pressio::rom::lspg::impl::composeUnsteady<
    ::pressio::rom::lspg::impl::PreconditionedHyperReduced, void, T1,
    typename std::remove_cv<typename std::remove_reference<Args>::type>::type...
    >,
  ::pressio::rom::lspg::impl::composeSteady<
    ::pressio::rom::lspg::impl::PreconditionedHyperReduced, void, T1,
    typename std::remove_cv<typename std::remove_reference<Args>::type>::type...
    >
  >::type;

template<typename T1, typename ...Args>
using composePreconditionedHyperReducedProblem_t =
  typename composePreconditionedHyperReducedProblem<T1, Args...>::type;


/*
  =====================
  === masked ===
  =====================

  steady:
  template<fom_type, decoder_t, romstate_t, masker_t>

  unsteady cont-time api:
  template<stepper_tag, fom_type, decoder_t, romstate_t, masker_t>

  unsteady discrete-time api:
  template<stepper_tag, fom_type, decoder_t, romstate_t, masker_t>
*/
template<typename T1, typename ...Args>
using composeMaskedProblem =
  typename std::conditional<
  ::pressio::ode::predicates::is_stepper_tag<T1>::value,
  ::pressio::rom::lspg::impl::composeUnsteady<
    ::pressio::rom::lspg::impl::Masked, void, T1,
    typename std::remove_cv<typename std::remove_reference<Args>::type>::type...
    >,
  ::pressio::rom::lspg::impl::composeSteady<
    ::pressio::rom::lspg::impl::Masked, void, T1,
    typename std::remove_cv<typename std::remove_reference<Args>::type>::type...
    >
  >::type;

template<typename T1, typename ...Args>
using composeMaskedProblem_t = typename composeMaskedProblem<T1, Args...>::type;

}}}}
#endif  // ROM_LSPG_IMPL_ROM_COMPOSE_PROBLEM_HPP_
