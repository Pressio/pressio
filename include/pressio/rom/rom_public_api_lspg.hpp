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

#ifndef ROM_LSPG_ROM_CREATE_DEFAULT_LSPG_PROBLEM_HPP_
#define ROM_LSPG_ROM_CREATE_DEFAULT_LSPG_PROBLEM_HPP_

#include "./impl/rom_lspg_steady_compose_impl.hpp"

namespace pressio{ namespace rom{ namespace lspg{

// --------------------
// for steady
// --------------------
template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ReturnType = typename impl::ComposeDefaultProblemSteady<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType
    >::type
  >
ReturnType create_default_steady_problem(const FomSystemType & fomSysObj,
					 DecoderType & decoder,
					 const RomStateType & stateIn,
					 const FomReferenceStateType & fomRef)
{
  return ReturnType(fomSysObj, decoder, stateIn, fomRef);
}

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class PreconditionerType,
  class ReturnType = typename impl::ComposePrecDefaultProblemSteady<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType, PreconditionerType
    >::type
  >
ReturnType create_default_steady_problem(const FomSystemType & fomSysObj,
					 DecoderType & decoder,
					 const RomStateType & stateIn,
					 const FomReferenceStateType & fomRef,
					 const PreconditionerType & prec)
{
  return ReturnType(fomSysObj, decoder, stateIn, fomRef, prec);
}

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class MaskerType,
  class ReturnType = typename impl::ComposeMaskedProblemSteady<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType, MaskerType
    >::type
  >
ReturnType create_masked_steady_problem(const FomSystemType & fomSysObj,
					DecoderType & decoder,
					const RomStateType & stateIn,
					const FomReferenceStateType & fomRef,
					const MaskerType & masker)
{
  return ReturnType(fomSysObj, decoder, stateIn, fomRef, masker);
}

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class MaskerType,
  class PreconditionerType,
  class ReturnType = typename impl::ComposePrecMaskedProblemSteady<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType,
    MaskerType, PreconditionerType
    >::type
  >
ReturnType create_masked_steady_problem(const FomSystemType & fomSysObj,
					DecoderType & decoder,
					const RomStateType & stateIn,
					const FomReferenceStateType & fomRef,
					const MaskerType & masker,
					const PreconditionerType & prec)
{
  return ReturnType(fomSysObj, decoder, stateIn, fomRef, prec, masker);
}

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class ReturnType = typename impl::ComposeHyperreducedProblemSteady<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType
    >::type
  >
ReturnType create_hyperreduced_steady_problem(const FomSystemType & fomSysObj,
					      DecoderType & decoder,
					      const RomStateType & stateIn,
					      const FomReferenceStateType & fomRef)
{
  return ReturnType(fomSysObj, decoder, stateIn, fomRef);
}

template<
  class FomSystemType,
  class DecoderType,
  class RomStateType,
  class FomReferenceStateType,
  class PreconditionerType,
  class ReturnType = typename impl::ComposePrecHypredProblemSteady<
    FomSystemType, DecoderType, RomStateType, FomReferenceStateType, PreconditionerType
    >::type
  >
ReturnType create_hyperreduced_steady_problem(const FomSystemType & fomSysObj,
					      DecoderType & decoder,
					      const RomStateType & stateIn,
					      const FomReferenceStateType & fomRef,
					      const PreconditionerType & prec)
{
  return ReturnType(fomSysObj, decoder, stateIn, fomRef, prec);
}

}}}//end namespace pressio::rom::lspg
#endif  // ROM_LSPG_ROM_CREATE_LSPG_PROBLEM_HPP_
