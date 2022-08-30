/*
//@HEADER
// ************************************************************************
//
// rom_public_api_galerkin.hpp
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

#ifndef ROM_ROM_PUBLIC_API_GALERKIN_HPP_
#define ROM_ROM_PUBLIC_API_GALERKIN_HPP_

#include "./impl/rom_galerkin_compose.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

// ============================
// for cont-time, explicit
// ============================
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ReturnType = impl::ComposeContTimeExplicit_t<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType
    >
  >
ReturnType create_default_explicit_problem(::pressio::ode::StepScheme name,
					   const FomSystemType & fomSysObj,
					   DecoderType & decoder,
					   const RomStateType & stateIn,
					   const FomReferenceStateType & fomRef)
{

  impl::ensure_explicit_or_throw("galerkin_default_explicit", name);
  return ReturnType(name, fomSysObj, decoder, stateIn, fomRef);
}

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ProjectorType,
  class ReturnType = impl::ComposeContTimeExplicit_t<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType, ProjectorType
    >
  >
ReturnType create_hyperreduced_explicit_problem(::pressio::ode::StepScheme name,
						const FomSystemType & fomSysObj,
						DecoderType & decoder,
						const RomStateType & stateIn,
						const FomReferenceStateType & fomRef,
						const ProjectorType & projector)
{

  impl::ensure_explicit_or_throw("galerkin_hyperreduced_explicit", name);
  return ReturnType(name, fomSysObj, decoder, stateIn, fomRef, projector);
}

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ProjectorType,
  class MaskerType,
  class ReturnType = impl::ComposeContTimeExplicit_t<
    FomSystemType, DecoderType, RomStateType,
    FomReferenceStateType, MaskerType, ProjectorType
    >
  >
ReturnType create_masked_explicit_problem(::pressio::ode::StepScheme name,
					  const FomSystemType & fomSysObj,
					  DecoderType & decoder,
					  const RomStateType & stateIn,
					  const FomReferenceStateType & fomRef,
					  const ProjectorType & projector,
					  const MaskerType & masker)
{

  impl::ensure_explicit_or_throw("galerkin_masked_explicit", name);
  return ReturnType(name, fomSysObj, decoder, stateIn, fomRef, projector, masker);
}

// ============================
// for cont-time, implicit
// ============================
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ReturnType = impl::ComposeContTimeImplicit_t<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType
    >
  >
ReturnType create_default_implicit_problem(::pressio::ode::StepScheme name,
					   const FomSystemType & fomSysObj,
					   DecoderType & decoder,
					   const RomStateType & stateIn,
					   const FomReferenceStateType & fomRef)
{

  impl::ensure_implicit_or_throw("galerkin_default_implicit", name);
  return ReturnType(name, fomSysObj, decoder, stateIn, fomRef);
}

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ProjectorType,
  class ReturnType = impl::ComposeContTimeImplicit_t<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType, ProjectorType
    >
  >
ReturnType create_hyperreduced_implicit_problem(::pressio::ode::StepScheme name,
						const FomSystemType & fomSysObj,
						DecoderType & decoder,
						const RomStateType & stateIn,
						const FomReferenceStateType & fomRef,
						const ProjectorType & projector)
{

  impl::ensure_implicit_or_throw("galerkin_hyperreduced_implicit", name);
  return ReturnType(name, fomSysObj, decoder, stateIn, fomRef, projector);
}

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ProjectorType,
  class MaskerType,
  class ReturnType = impl::ComposeContTimeImplicit_t<
    FomSystemType, DecoderType, RomStateType,
    FomReferenceStateType, MaskerType, ProjectorType
    >
  >
ReturnType create_masked_implicit_problem(::pressio::ode::StepScheme name,
					  const FomSystemType & fomSysObj,
					  DecoderType & decoder,
					  const RomStateType & stateIn,
					  const FomReferenceStateType & fomRef,
					  const ProjectorType & projector,
					  const MaskerType & masker)
{

  impl::ensure_implicit_or_throw("galerkin_masked_implicit", name);
  return ReturnType(name, fomSysObj, decoder, stateIn, fomRef, projector, masker);
}

// ============================
// for discrete-time
// ============================
template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ReturnType = impl::ComposeDiscTime_t<
    num_states, FomSystemType, DecoderType, RomStateType, FomReferenceStateType
    >
  >
ReturnType create_default_problem(const FomSystemType & fomSysObj,
				  DecoderType & decoder,
				  const RomStateType & stateIn,
				  const FomReferenceStateType & fomRef)
{
  return ReturnType(::pressio::ode::StepScheme::ImplicitArbitrary,
		    fomSysObj, decoder, stateIn, fomRef);
}

template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ProjectorType,
  class ReturnType = impl::ComposeDiscTime_t<
    num_states, FomSystemType, DecoderType, RomStateType, FomReferenceStateType, ProjectorType
    >
  >
ReturnType create_hyperreduced_problem(const FomSystemType & fomSysObj,
				       DecoderType & decoder,
				       const RomStateType & stateIn,
				       const FomReferenceStateType & fomRef,
				       const ProjectorType & projector)
{

  return ReturnType(::pressio::ode::StepScheme::ImplicitArbitrary,
		    fomSysObj, decoder, stateIn, fomRef, projector);
}

template<
  std::size_t num_states,
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ProjectorType,
  class MaskerType,
  class ReturnType = impl::ComposeDiscTime_t<
    num_states, FomSystemType, DecoderType,
    RomStateType, FomReferenceStateType, ProjectorType, MaskerType
    >
  >
ReturnType create_masked_problem(const FomSystemType & fomSysObj,
				 DecoderType & decoder,
				 const RomStateType & stateIn,
				 const FomReferenceStateType & fomRef,
				 const ProjectorType & projector,
				 const MaskerType & masker)
{

  return ReturnType(::pressio::ode::StepScheme::ImplicitArbitrary,
		    fomSysObj, decoder, stateIn, fomRef, projector, masker);
}

}}}//end namespace pressio::rom::galerkin
#endif  // ROM_ROM_PUBLIC_API_GALERKIN_HPP_
