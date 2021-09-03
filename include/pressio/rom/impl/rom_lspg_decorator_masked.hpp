/*
//@HEADER
// ************************************************************************
//
// rom_masked.hpp
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

#ifndef ROM_LSPG_DECORATORS_ROM_MASKED_HPP_
#define ROM_LSPG_DECORATORS_ROM_MASKED_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template <class DataType, class MaskerType, class MaskableType>
class Masked : public MaskableType
{

public:
  Masked() = delete;
  Masked(const Masked &) = default;
  Masked & operator=(const Masked &) = default;
  Masked(Masked &&) = default;
  Masked & operator=(Masked &&) = default;
  ~Masked() = default;

  template <class ... Args>
  Masked(const MaskerType & maskerObj, Args && ... args)
    : MaskableType(std::forward<Args>(args)...),
      unmaskedObject_(MaskableType::create()),
      masker_(maskerObj)
  {}

public:
  DataType create() const{
    return DataType(masker_.get().createApplyMaskResult(unmaskedObject_));
  }

  // steady calls this
  template <class LspgStateType>
  void compute(const LspgStateType & state, DataType & maskedResult) const
  {
    MaskableType::compute(state, unmaskedObject_);
    masker_(unmaskedObject_, maskedResult);
  }

  // // unsteady calls this
  // template <
  //   class stepper_tag,
  //   class LspgStateType,
  //   class prev_states_t,
  //   class fom_system_t,
  //   class time_type
  //   >
  // void compute(const LspgStateType & state,
  // 	       const prev_states_t & prevStates,
  // 	       const fom_system_t & systemObj,
  // 	       const time_type & time,
  // 	       const time_type & dt,
  // 	       const ::pressio::ode::step_count_type & step,
  // 	       DataType & maskedResult) const
  // {
  //   MaskableType::template compute<stepper_tag>
  //     (state, prevStates, systemObj, time, dt, step, unmaskedObject_);
  //   masker_.get().applyMask(*unmaskedObject_.data(), time, *maskedResult.data());
  // }

private:
  mutable DataType unmaskedObject_;
  std::reference_wrapper<const MaskerType> masker_;

};

}}}}
#endif  // ROM_LSPG_DECORATORS_ROM_MASKED_HPP_
