/*
//@HEADER
// ************************************************************************
//
// rom_create_default_galerkin_problem.hpp
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

#ifndef ROM_GALERKIN_ROM_CREATE_DEFAULT_GALERKIN_PROBLEM_HPP_
#define ROM_GALERKIN_ROM_CREATE_DEFAULT_GALERKIN_PROBLEM_HPP_

#include "./impl/rom_galerkin_cont_time_compose_impl.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

template<
  class StepperTag,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ReturnType = impl::ComposeDefaultProblemContTime_t<
    StepperTag, void, FomSystemType, DecoderType, RomStateType, FomReferenceStateType
    >
  >
mpl::enable_if_t<::pressio::ode::is_stepper_tag<StepperTag>::value, ReturnType>
create_default_problem(const FomSystemType & fomSysObj,
                       DecoderType & decoder,
                       const RomStateType & stateIn,
                       const FomReferenceStateType & fomRef)
{
  return ReturnType(fomSysObj, decoder, stateIn, fomRef);
}

template<
  class StepperTag,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ProjectorType,
  class MaskerType,
  class ReturnType = impl::ComposeMaskedVelocityProblemContTime_t<
    StepperTag, void, FomSystemType, DecoderType, RomStateType,
    FomReferenceStateType, MaskerType, ProjectorType
    >
  >
mpl::enable_if_t<::pressio::ode::is_stepper_tag<StepperTag>::value, ReturnType>
create_masked_velocity_problem(const FomSystemType & fomSysObj,
                               DecoderType & decoder,
                               const RomStateType & stateIn,
                               const FomReferenceStateType & fomRef,
                               const ProjectorType & projector,
                               const MaskerType & masker)

{
  return ReturnType(fomSysObj, decoder, stateIn, fomRef, projector, masker);
}

template<
  class StepperTag,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ProjectorType,
  class ReturnType = impl::ComposeHypRedVeloProblemContTime_t<
    StepperTag, void, FomSystemType, DecoderType, RomStateType,
    FomReferenceStateType, ProjectorType
    >
  >
mpl::enable_if_t<::pressio::ode::is_stepper_tag<StepperTag>::value, ReturnType>
create_hyperreduced_velocity_problem(const FomSystemType & fomSysObj,
				     DecoderType & decoder,
				     const RomStateType & stateIn,
				     const FomReferenceStateType & fomRef,
				     const ProjectorType & projector)

{
  return ReturnType(fomSysObj, decoder, stateIn, fomRef, projector);
}

}}}//end namespace pressio::rom::galerkin
#endif  // ROM_GALERKIN_ROM_CREATE_DEFAULT_GALERKIN_PROBLEM_HPP_




// // Note that here the user specifies the galerkin jacobian type
// template<
//   typename rom_jacobian_type,
//   std::size_t order,
//   std::size_t totNumStates,
//   typename fom_system_type,
//   typename decoder_type,
//   typename rom_state_type,
//   typename fom_native_state,
//   typename return_t =
//     impl::composeDefaultProblemDiscTime_t<
//       ::pressio::ode::implicitmethods::Arbitrary,
//       fom_system_type, decoder_type, rom_state_type, rom_jacobian_type,
//       ::pressio::ode::types::StepperOrder<order>,
//       ::pressio::ode::types::StepperTotalNumberOfStates<totNumStates>
//       >
//   >
// return_t
// createDefaultProblem(const fom_system_type & fomSysObj,
//          decoder_type & decoder,
//          const rom_state_type & stateIn,
//          const fom_native_state & fomRef)
// {
//   return return_t(fomSysObj, decoder, stateIn, fomRef);
// }

// // Note that here the user does NOT specify the galerkin jacobian type
// // so pressio has some rules behind the scenes to select the jacobian type
// // based on the galerkin state
// template<
//   std::size_t order,
//   std::size_t totNumStates,
//   typename fom_system_type,
//   typename decoder_type,
//   typename rom_state_type,
//   typename fom_native_state,
//   typename return_t =
//     impl::composeDefaultProblemDiscTime_t<
//       ::pressio::ode::implicitmethods::Arbitrary,
//       fom_system_type, decoder_type, rom_state_type, void,
//       ::pressio::ode::types::StepperOrder<order>,
//       ::pressio::ode::types::StepperTotalNumberOfStates<totNumStates>
//       >
//   >
// mpl::enable_if_t<
//   ::pressio::rom::constraints::most_likely_discrete_time_system<fom_system_type>::value,
//   return_t
//   >
// createDefaultProblem(const fom_system_type & fomSysObj,
//          decoder_type & decoder,
//          const rom_state_type & stateIn,
//          const fom_native_state & fomRef)
// {
//   return return_t(fomSysObj, decoder, stateIn, fomRef);
// }




// OLD API (before Miko edits)
// template<
//   typename stepper_tag,
//   typename fom_system_type,
//   typename decoder_type,
//   typename rom_state_type,
//   typename fom_native_state,
//   typename ...Args
//   >
// mpl::enable_if_t<
//   ::pressio::rom::constraints::most_likely_continuous_time_system<fom_system_type>::value,
//   impl::composeDefaultProblemContTime_t<stepper_tag, fom_system_type, decoder_type, rom_state_type, Args...>
//   >
// createDefaultProblem(const fom_system_type & fomSysObj,
//         decoder_type & decoder,
//         const rom_state_type & stateIn,
//         const fom_native_state & fomRef,
//         Args && ... args)
// {
//   using return_t = impl::composeDefaultProblemContTime_t<
//     stepper_tag, fom_system_type, decoder_type, rom_state_type, Args...>;

//   return return_t(fomSysObj, decoder, stateIn,
//      fomRef, std::forward<Args>(args)...);
// }

// // Note that here the user specifies the galerkin jacobian type
// template<
//   typename rom_jacobian_type,
//   std::size_t order,
//   std::size_t totNumStates,
//   typename fom_system_type,
//   typename decoder_type,
//   typename rom_state_type,
//   typename fom_native_state
//   >
// mpl::enable_if_t<
//   ::pressio::rom::constraints::most_likely_discrete_time_system<fom_system_type>::value,
//   impl::composeDefaultProblemDiscTime_t<
//     pressio::ode::ImplicitArbitrary,
//     fom_system_type, decoder_type, rom_state_type, rom_jacobian_type,
//     ::pressio::ode::StepperOrder<order>,
//     ::pressio::ode::StepperTotalNumberOfStates<totNumStates>
//     >
//   >
// createDefaultProblem(const fom_system_type & fomSysObj,
//         decoder_type & decoder,
//         const rom_state_type & stateIn,
//         const fom_native_state & fomRef)
// {
//   using return_t =
//     impl::composeDefaultProblemDiscTime_t<
//       ::pressio::ode::ImplicitArbitrary,
//     fom_system_type, decoder_type, rom_state_type, rom_jacobian_type,
//     ::pressio::ode::StepperOrder<order>,
//     ::pressio::ode::StepperTotalNumberOfStates<totNumStates>
//     >;
//   return return_t(fomSysObj, decoder, stateIn, fomRef);
// }

// // Note that here the user does NOT specify the galerkin jacobian type
// // so pressio has some rules behind the scenes to select the jacobian type
// // based on the galerkin state
// template<
//   std::size_t order,
//   std::size_t totNumStates,
//   typename fom_system_type,
//   typename decoder_type,
//   typename rom_state_type,
//   typename fom_native_state
//   >
// mpl::enable_if_t<
//   ::pressio::rom::constraints::most_likely_discrete_time_system<fom_system_type>::value,
//   impl::composeDefaultProblemDiscTime_t<
//     pressio::ode::ImplicitArbitrary,
//     fom_system_type, decoder_type, rom_state_type, void,
//     ::pressio::ode::StepperOrder<order>,
//     ::pressio::ode::StepperTotalNumberOfStates<totNumStates>
//     >
//   >
// createDefaultProblem(const fom_system_type & fomSysObj,
//         decoder_type & decoder,
//         const rom_state_type & stateIn,
//         const fom_native_state & fomRef)
// {
//   using return_t =
//     impl::composeDefaultProblemDiscTime_t<
//       ::pressio::ode::ImplicitArbitrary,
//     fom_system_type, decoder_type, rom_state_type, void,
//     ::pressio::ode::StepperOrder<order>,
//     ::pressio::ode::StepperTotalNumberOfStates<totNumStates>
//     >;
//   return return_t(fomSysObj, decoder, stateIn, fomRef);
// }
