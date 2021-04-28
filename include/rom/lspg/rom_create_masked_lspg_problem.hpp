/*
//@HEADER
// ************************************************************************
//
// rom_create_masked_lspg_problem.hpp
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

#ifndef ROM_LSPG_ROM_CREATE_MASKED_LSPG_PROBLEM_HPP_
#define ROM_LSPG_ROM_CREATE_MASKED_LSPG_PROBLEM_HPP_

namespace pressio{ namespace rom{ namespace lspg{

// create masked steady
template<
  typename fom_system_type,
  typename decoder_type,
  typename rom_state_type,
  typename fom_native_state,
  typename masker_type
  >
mpl::enable_if_t<
  ::pressio::rom::constraints::most_likely_steady_system<fom_system_type>::value,
  impl::composeMaskedProblem_t<
    fom_system_type, decoder_type, rom_state_type, masker_type
    >
  >
createMaskedProblemSteady(const fom_system_type & fomSysObj,
			  decoder_type & decoder,
			  const rom_state_type & romStateIn,
			  const fom_native_state & fomRef,
			  const masker_type & masker)
{
  using return_t = impl::composeMaskedProblem_t<
    fom_system_type, decoder_type, rom_state_type, masker_type>;

  static_assert
    (std::is_same<fom_native_state, typename return_t::fom_native_state_t>::value,
     "The fom reference state type deduced for the create function is not \
compatible with the fom state type detected from adapter class");

  return return_t(fomSysObj, decoder, romStateIn,
		  fomRef, masker);
}

// unsteady (continuous-time api)
template<
  typename odetag,
  typename fom_system_type,
  typename decoder_type,
  typename rom_state_type,
  typename fom_native_state,
  typename masker_type
  >
mpl::enable_if_t<
  ::pressio::rom::constraints::most_likely_continuous_time_system<fom_system_type>::value,
  impl::composeMaskedProblem_t<
    odetag, fom_system_type, decoder_type, rom_state_type, masker_type
    >
  >
createMaskedProblemUnsteady(const fom_system_type & fomSysObj,
			    decoder_type & decoder,
			    const rom_state_type & romStateIn,
			    const fom_native_state & fomRef,
			    const masker_type & masker)
{
  using return_t = impl::composeMaskedProblem_t<
    odetag, fom_system_type, decoder_type, rom_state_type, masker_type>;

  static_assert
    (std::is_same<fom_native_state, typename return_t::fom_native_state_t>::value,
     "The fom reference state type deduced for the create function is not \
compatible with the fom state type detected from adapter class");

  return return_t(fomSysObj, decoder, romStateIn, fomRef, masker);
}

//  unsteady (discrete-time api)
template<
  std::size_t order,
  std::size_t totNumStates,
  typename fom_system_type,
  typename decoder_type,
  typename rom_state_type,
  typename fom_native_state,
  typename masker_type
  >
mpl::enable_if_t<
  ::pressio::rom::constraints::most_likely_discrete_time_system<fom_system_type>::value,
  impl::composeMaskedProblem_t<
    pressio::ode::implicitmethods::Arbitrary,
    fom_system_type, decoder_type, rom_state_type, masker_type,
    ::pressio::ode::types::StepperOrder<order>,
    ::pressio::ode::types::StepperTotalNumberOfStates<totNumStates>
    >
  >
createMaskedProblemUnsteady(const fom_system_type & fomSysObj,
			    decoder_type & decoder,
			    const rom_state_type & romStateIn,
			    const fom_native_state & fomNominalState,
			    const masker_type & masker)
{
  using return_t = impl::composeMaskedProblem_t<
    pressio::ode::implicitmethods::Arbitrary,
    fom_system_type, decoder_type, rom_state_type, masker_type,
    ::pressio::ode::types::StepperOrder<order>,
    ::pressio::ode::types::StepperTotalNumberOfStates<totNumStates>
    >;

  static_assert
    (std::is_same<fom_native_state, typename return_t::fom_native_state_t>::value,
     "The type deduced for the FOM nominal state passed to the create function is not \
compatible with the FOM state type detected from adapter class");

  return return_t(fomSysObj, decoder, romStateIn, fomNominalState, masker);
}

}}}
#endif  // ROM_LSPG_ROM_CREATE_MASKED_LSPG_PROBLEM_HPP_
