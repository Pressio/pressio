/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_problem_generator.hpp
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

#ifndef ROM_GALERKIN_PROBLEM_GENERATOR_HPP_
#define ROM_GALERKIN_PROBLEM_GENERATOR_HPP_

#include "rom_galerkin_problem_type_generator_default.hpp"

namespace pressio{ namespace rom{

template <typename problem_t>
class GalerkinProblemGenerator<problem_t>
  : public problem_t
{
public:
  using typename problem_t::fom_t;
  using typename problem_t::scalar_t;
  using typename problem_t::fom_native_state_t;
  using typename problem_t::fom_state_t;
  using typename problem_t::fom_velocity_t;

  using typename problem_t::galerkin_state_t;
  using typename problem_t::decoder_t;
  using typename problem_t::fom_state_reconstr_t;
  using typename problem_t::fom_states_data;
  using typename problem_t::ud_ops_t;

  using typename problem_t::galerkin_residual_policy_t;
  using typename problem_t::galerkin_stepper_t;

private:
  fom_state_t			fomStateReference_;
  fom_state_reconstr_t		fomStateReconstructor_;
  fom_velocity_t		fomVelocityRef_;
  fom_states_data		fomStates_;
  galerkin_residual_policy_t	residualPolicy_;
  galerkin_stepper_t		stepperObj_;

public:
  galerkin_stepper_t & getStepperRef(){
    return stepperObj_;
  }

  const fom_state_reconstr_t & getFomStateReconstructorCRef() const{
    return fomStateReconstructor_;
  }

public:
  /*
   * ud_ops_t == void and state_type is a wrapper
  */
  template <
    typename _ud_ops_t = ud_ops_t,
    typename ::pressio::mpl::enable_if_t<
      std::is_void<_ud_ops_t>::value and
      ::pressio::containers::meta::is_wrapper<galerkin_state_t>::value
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
      and !::pressio::containers::meta::is_array_pybind11<galerkin_state_t>::value
#endif
      > * = nullptr
  >
  GalerkinProblemGenerator(const fom_t		    & appObj,
  			   const fom_native_state_t & yFomRefNative,
  			   const decoder_t	    & decoder,
  			   galerkin_state_t	    & yROM,
  			   scalar_t		    t0)
    : fomStateReference_(yFomRefNative),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomVelocityRef_( appObj.velocity(*fomStateReference_.data(), t0) ),
      fomStates_(fomStateReconstructor_, fomStateReference_),
      residualPolicy_(fomVelocityRef_, fomStates_, decoder),
      stepperObj_(yROM, appObj, residualPolicy_)
  {}

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  /*
   * ud_ops_t == void and state_type is pybind11::array
  */
  template <
    typename _ud_ops_t = ud_ops_t,
    ::pressio::mpl::enable_if_t<
      std::is_void<_ud_ops_t>::value and
      ::pressio::containers::meta::is_array_pybind11<galerkin_state_t>::value
      > * = nullptr
  >
  GalerkinProblemGenerator(const fom_t		    & appObj,
  			   const fom_native_state_t & yFomRefNative,
  			   const decoder_t	    & decoder,
  			   galerkin_state_t	    & yROM,
  			   scalar_t		    t0)
    : fomStateReference_(yFomRefNative),
      fomStateReconstructor_(fomStateReference_, decoder),
      fomVelocityRef_( appObj.attr("velocity")(fomStateReference_, t0) ),
      fomStates_(fomStateReconstructor_, fomStateReference_),
      residualPolicy_(fomVelocityRef_, fomStates_, decoder),
      stepperObj_(yROM, appObj, residualPolicy_)
  {}
#endif

};

}}//end namespace pressio::rom
#endif
