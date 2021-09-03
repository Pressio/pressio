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

template <class DataType, class PreconditionerType, class preconditionable_policy>
class Preconditioned : public preconditionable_policy
{
public:
  Preconditioned() = delete;
  Preconditioned(const Preconditioned &) = default;
  Preconditioned & operator=(const Preconditioned &) = default;
  Preconditioned(Preconditioned &&) = default;
  Preconditioned & operator=(Preconditioned &&) = default;
  ~Preconditioned() = default;

  Preconditioned(const preconditionable_policy & obj)
    : preconditionable_policy(obj)
  {}

  template <class ... Args>
  Preconditioned(const PreconditionerType & preconditionerIn,
		 Args && ... args)
    : preconditionable_policy(std::forward<Args>(args)...),
      preconditionerObj_(preconditionerIn)
  {}

public:
  DataType create() const{ return preconditionable_policy::create(); }

  // steady calls this
  template <class LspgStateType>
  void compute(const LspgStateType & currentState, DataType & unpreconditionedObj) const
  {
    preconditionable_policy::compute(currentState, unpreconditionedObj);
    const auto & fomState = fomStatesMngr_.get().currentFomState();
    preconditionerObj_(fomState, unpreconditionedObj);
  }

  // // unsteady calls this
  // template <
  //   class stepper_tag,
  //   class prev_states_t,
  //   class fom_system_t,
  //   class scalar_t,
  //   class state_t
  //   >
  // void compute(const state_t & currentState,
  // 	       const prev_states_t & prevStates,
  // 	       const fom_system_t & systemObj,
  // 	       const scalar_t & time,
  // 	       const scalar_t & dt,
  // 	       const ::pressio::ode::step_count_type & step,
  // 	       DataType & unpreconditionedField) const
  // {
  //   preconditionable_policy::template compute<
  //     stepper_tag>(currentState, prevStates, systemObj,
  //       time, dt, step, unpreconditionedField);

  //   const auto & yFom = fomStatesMngr_(::pressio::ode::nPlusOne());
  //   preconditionerObj_.get().applyPreconditioner(*yFom.data(), time,
  // 						 *unpreconditionedField.data());
  // }

private:
  using preconditionable_policy::fomStatesMngr_;
  std::reference_wrapper<const PreconditionerType> preconditionerObj_;

};

}}}} //end namespace pressio::rom::lspg::decorator
#endif  // ROM_LSPG_DECORATORS_ROM_PRECONDITIONED_HPP_
