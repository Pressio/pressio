/*
//@HEADER
// ************************************************************************
//
// rom_compose_galerkin.hpp
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

#ifndef ROM_GALERKIN_ROM_COMPOSE_GALERKIN_HPP_
#define ROM_GALERKIN_ROM_COMPOSE_GALERKIN_HPP_

#include "./impl/rom_compose_galerkin_impl.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

template<typename ...Args>
using composeDefaultProblem = impl::compose<
  impl::Default, void,
  typename std::remove_cv<typename std::remove_reference<Args>::type>::type...>;

template<typename ...Args>
using composeDefaultProblem_t = typename composeDefaultProblem<Args...>::type;

template<
  typename tag,
  typename fom_system_type,
  typename decoder_type,
  typename rom_state_type,
  typename fom_native_state,
  typename ...Args
  >
mpl::enable_if_t<
  ::pressio::rom::concepts::continuous_time_explicit_system<fom_system_type>::value,
  galerkin::composeDefaultProblem_t<
    tag, fom_system_type, decoder_type, rom_state_type, Args...
    >
  >
createDefaultProblem(const fom_system_type & fomSysObj,
		     const decoder_type & decoder,
		     const rom_state_type & stateIn,
		     const fom_native_state & fomRef,
		     Args && ... args)
{
  using return_t = galerkin::composeDefaultProblem_t<
    tag, fom_system_type, decoder_type, rom_state_type, Args...>;

  static_assert
  (std::is_same<fom_native_state, typename return_t::fom_native_state_t>::value,
   "The fom reference state type deduced for the create function is not \
compatible with the fom state type detected from adapter class");

  return return_t(fomSysObj, decoder, stateIn,
		  fomRef, std::forward<Args>(args)...);
}


template<
  typename rom_jacobian_type,
  std::size_t order,
  std::size_t totNumStates,
  typename fom_system_type,
  typename decoder_type,
  typename rom_state_type,
  typename fom_native_state
  >
mpl::enable_if_t<
  ::pressio::rom::concepts::discrete_time_system<fom_system_type>::value,
  galerkin::composeDefaultProblem_t<
    pressio::ode::implicitmethods::Arbitrary,
    fom_system_type, decoder_type, rom_state_type, rom_jacobian_type,
    ::pressio::ode::types::StepperOrder<order>,
    ::pressio::ode::types::StepperTotalNumberOfStates<totNumStates>
    >
  >
createDefaultProblem(const fom_system_type & fomSysObj,
		     const decoder_type & decoder,
		     const rom_state_type & stateIn,
		     const fom_native_state & fomRef)
{
  using return_t =
    galerkin::composeDefaultProblem_t<
      ::pressio::ode::implicitmethods::Arbitrary,
    fom_system_type, decoder_type, rom_state_type, rom_jacobian_type,
    ::pressio::ode::types::StepperOrder<order>,
    ::pressio::ode::types::StepperTotalNumberOfStates<totNumStates>
    >;

  static_assert
  (std::is_same<fom_native_state, typename return_t::fom_native_state_t>::value,
   "The fom reference state type deduced for the create function is not \
compatible with the fom state type detected from adapter class");

  return return_t(fomSysObj, decoder, stateIn, fomRef);
}

}}}//end namespace pressio::rom::galerkin
#endif  // ROM_GALERKIN_ROM_COMPOSE_GALERKIN_HPP_
