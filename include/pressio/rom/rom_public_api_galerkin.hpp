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

#include "./impl/rom_galerkin_compose_impl.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

//
// for cont-time
//
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
  class ReturnType = impl::ComposeHypRedVeloProblemContTime_t<
    StepperTag, void, FomSystemType, DecoderType, RomStateType,
    FomReferenceStateType, ProjectorType
    >
  >
mpl::enable_if_t<::pressio::ode::is_stepper_tag<StepperTag>::value, ReturnType>
create_hyperreduced_problem(const FomSystemType & fomSysObj,
			    DecoderType & decoder,
			    const RomStateType & stateIn,
			    const FomReferenceStateType & fomRef,
			    const ProjectorType & projector)

{
  return ReturnType(fomSysObj, decoder, stateIn, fomRef, projector);
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
create_masked_problem(const FomSystemType & fomSysObj,
		      DecoderType & decoder,
		      const RomStateType & stateIn,
		      const FomReferenceStateType & fomRef,
		      const ProjectorType & projector,
		      const MaskerType & masker)

{
  return ReturnType(fomSysObj, decoder, stateIn, fomRef, projector, masker);
}

//
// for discrete-time
//
template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ReturnType = typename impl::ComposeDefaultProblemDiscTime<
    num_states, FomSystemType, DecoderType, RomStateType, FomReferenceStateType
    >::type
  >
ReturnType create_default_problem(const FomSystemType & fomSysObj,
				  DecoderType & decoder,
				  const RomStateType & stateIn,
				  const FomReferenceStateType & fomRef)
{
  return ReturnType(fomSysObj, decoder, stateIn, fomRef);
}

template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ProjectorType,
  class ReturnType = typename impl::ComposeHypRedProblemDiscTime<
    num_states, FomSystemType, DecoderType, RomStateType, FomReferenceStateType, ProjectorType
    >::type
  >
ReturnType create_hyperreduced_problem(const FomSystemType & fomSysObj,
				       DecoderType & decoder,
				       const RomStateType & stateIn,
				       const FomReferenceStateType & fomRef,
				       const ProjectorType & projector)

{
  return ReturnType(fomSysObj, decoder, stateIn, fomRef, projector);
}

template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ProjectorType,
  class MaskerType,
  class ReturnType = typename impl::ComposeMaskedProblemDiscTime<
    num_states, FomSystemType, DecoderType,
    RomStateType, FomReferenceStateType, ProjectorType, MaskerType
    >::type
  >
ReturnType create_masked_problem(const FomSystemType & fomSysObj,
				 DecoderType & decoder,
				 const RomStateType & stateIn,
				 const FomReferenceStateType & fomRef,
				 const ProjectorType & projector,
				 const MaskerType & masker)

{
  return ReturnType(fomSysObj, decoder, stateIn, fomRef, projector, masker);
}

}}}//end namespace pressio::rom::galerkin
#endif  // ROM_GALERKIN_ROM_CREATE_DEFAULT_GALERKIN_PROBLEM_HPP_
