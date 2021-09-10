/*
//@HEADER
// ************************************************************************
//
// rom_preconditioned.hpp
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

#ifndef ROM_LSPG_DECORATORS_ROM_PRECONDITIONED_HPP_
#define ROM_LSPG_DECORATORS_ROM_PRECONDITIONED_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template <class DataType, class PreconditionerType, class preconditionable>
class PreconditionDecoratorSteady : public preconditionable
{
public:
  PreconditionDecoratorSteady() = delete;
  PreconditionDecoratorSteady(const PreconditionDecoratorSteady &) = default;
  PreconditionDecoratorSteady & operator=(const PreconditionDecoratorSteady &) = default;
  PreconditionDecoratorSteady(PreconditionDecoratorSteady &&) = default;
  PreconditionDecoratorSteady & operator=(PreconditionDecoratorSteady &&) = default;
  ~PreconditionDecoratorSteady() = default;

  PreconditionDecoratorSteady(const preconditionable & obj)
    : preconditionable(obj)
  {}

  template <class ... Args>
  PreconditionDecoratorSteady(const PreconditionerType & preconditionerIn,
			      Args && ... args)
    : preconditionable(std::forward<Args>(args)...),
      preconditionerObj_(preconditionerIn)
  {}

public:
  DataType create() const{ return preconditionable::create(); }

  template <class LspgStateType>
  void compute(const LspgStateType & currentState,
	       DataType & unpreconditionedObj) const
  {
    preconditionable::compute(currentState, unpreconditionedObj);
    const auto & fomState = fomStatesMngr_.get().currentFomState();
    preconditionerObj_(fomState, unpreconditionedObj);
  }

private:
  using preconditionable::fomStatesMngr_;
  std::reference_wrapper<const PreconditionerType> preconditionerObj_;
};


template <class DataType, class MaskerType, class MaskableType>
class MaskDecoratorSteady : public MaskableType
{

public:
  MaskDecoratorSteady() = delete;
  MaskDecoratorSteady(const MaskDecoratorSteady &) = default;
  MaskDecoratorSteady & operator=(const MaskDecoratorSteady &) = default;
  MaskDecoratorSteady(MaskDecoratorSteady &&) = default;
  MaskDecoratorSteady & operator=(MaskDecoratorSteady &&) = default;
  ~MaskDecoratorSteady() = default;

  template <class ... Args>
  MaskDecoratorSteady(const MaskerType & maskerObj, Args && ... args)
    : MaskableType(std::forward<Args>(args)...),
      unmaskedObject_(MaskableType::create()),
      masker_(maskerObj)
  {}

public:
  DataType create() const{
    return DataType(masker_.get().createApplyMaskResult(unmaskedObject_));
  }

  template <class LspgStateType>
  void compute(const LspgStateType & state, DataType & maskedResult) const
  {
    MaskableType::compute(state, unmaskedObject_);
    masker_(unmaskedObject_, maskedResult);
  }

private:
  mutable DataType unmaskedObject_;
  std::reference_wrapper<const MaskerType> masker_;
};

}}}} //end namespace pressio::rom::lspg::decorator
#endif  // ROM_LSPG_DECORATORS_ROM_PRECONDITIONED_HPP_
