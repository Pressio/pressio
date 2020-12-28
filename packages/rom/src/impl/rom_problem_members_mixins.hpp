/*
//@HEADER
// ************************************************************************
//
// rom_problem_members_mixins.hpp
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

#ifndef ROM_IMPL_ROM_PROBLEM_MEMBERS_MIXINS_HPP_
#define ROM_IMPL_ROM_PROBLEM_MEMBERS_MIXINS_HPP_

namespace pressio{ namespace rom{ namespace impl{

template<typename fom_system_t, bool isbinding=false>
struct FomObjMixin;

template<typename fom_system_t>
struct FomObjMixin<fom_system_t, false>
{
  std::reference_wrapper<const fom_system_t> fomObj_;

  FomObjMixin() = delete;
  FomObjMixin(const FomObjMixin &) = default;
  FomObjMixin & operator=(const FomObjMixin &) = delete;
  FomObjMixin(FomObjMixin &&) = default;
  FomObjMixin & operator=(FomObjMixin &&) = delete;
  ~FomObjMixin() = default;

  explicit FomObjMixin(const fom_system_t & fomObjIn)
    : fomObj_(fomObjIn){}

  const fom_system_t & fomCRef() const{ return fomObj_.get(); }
};

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template<typename fom_system_t>
struct FomObjMixin<fom_system_t, true>
{
  // when dealing with bindings for pressio4py, the fom_system_t
  // is a C++ class in pressio4py that wraps the actual FOM python object.
  // to construct this ROM problem, the Python code passes the
  // python FOM object NOT a C++ object instantiated from fom_system_t.
  // Therefore, ONLY when we deal with pressio4py, we create the fom obj
  // instead of referencing it.
  fom_system_t fomObj_;

  FomObjMixin() = delete;
  FomObjMixin(const FomObjMixin &) = default;
  FomObjMixin & operator=(const FomObjMixin &) = delete;
  FomObjMixin(FomObjMixin &&) = default;
  FomObjMixin & operator=(FomObjMixin &&) = delete;
  ~FomObjMixin() = default;

  explicit FomObjMixin(pybind11::object pyFomObj)
    : fomObj_(pyFomObj){}

  const fom_system_t & fomCRef() const{ return fomObj_; }
};
#endif

//---------------------------------------------------
template <
  typename T,
  typename ops_t,
  typename fom_state_t,
  typename fom_state_reconstr_t,
  typename fom_states_manager_t
  >
struct FomStatesMngrMixin : T
{
  const fom_state_t	     fomNominalState_;
  const fom_state_reconstr_t fomStateReconstructor_;
  fom_states_manager_t	     fomStatesMngr_;

  FomStatesMngrMixin() = delete;
  FomStatesMngrMixin(const FomStatesMngrMixin &) = default;
  FomStatesMngrMixin & operator=(const FomStatesMngrMixin &) = delete;
  FomStatesMngrMixin(FomStatesMngrMixin &&) = default;
  FomStatesMngrMixin & operator=(FomStatesMngrMixin &&) = delete;
  ~FomStatesMngrMixin() = default;

  template<
    typename T1, typename T2, typename T3, typename T4,
    typename _ops_t = ops_t,
    mpl::enable_if_t<std::is_void<_ops_t>::value, int > = 0
    >
  FomStatesMngrMixin(const T1 & fomObj,
		     const T2 & decoder,
		     const T3 & romStateIn,
		     const T4 & fomNominalStateNative)
    : T(fomObj),
      fomNominalState_(fomNominalStateNative),
      fomStateReconstructor_(fomNominalState_, decoder),
      fomStatesMngr_(fomStateReconstructor_, fomNominalState_)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(romStateIn);
  }

  template<
    typename T1, typename T2, typename T3, typename T4,
    typename _ops_t = ops_t,
    mpl::enable_if_t<!std::is_void<_ops_t>::value, int > = 0
    >
  FomStatesMngrMixin(const T1 & fomObj,
		     const T2 & decoder,
		     const T3 & romStateIn,
		     const T4 & fomNominalStateNative,
		     const _ops_t & ops)
    : T(fomObj),
      fomNominalState_(fomNominalStateNative),
      fomStateReconstructor_(fomNominalState_, decoder, ops),
      fomStatesMngr_(fomStateReconstructor_, fomNominalState_, ops)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(romStateIn);
  }
};

//---------------------------------------------------
// implicit stepper mixin
//---------------------------------------------------
// aux_stepper_t is valid
template <typename T, typename aux_stepper_t, typename stepper_t>
struct ImplicitStepperMixin : T
{
  aux_stepper_t auxStepperObj_;
  stepper_t stepperObj_;

  ImplicitStepperMixin() = delete;
  ImplicitStepperMixin(const ImplicitStepperMixin &) = default;
  ImplicitStepperMixin & operator=(const ImplicitStepperMixin &) = delete;
  ImplicitStepperMixin(ImplicitStepperMixin &&) = default;
  ImplicitStepperMixin & operator=(ImplicitStepperMixin &&) = delete;
  ~ImplicitStepperMixin() = default;

  template<typename T1, typename...Args>
  ImplicitStepperMixin(const T1 & romStateIn,
		       Args && ...args)
    : T(romStateIn, std::forward<Args>(args)...),
      auxStepperObj_(romStateIn, T::fomCRef(),
		     T::residualPolicy_, T::jacobianPolicy_),
      stepperObj_(romStateIn, T::fomCRef(),
		  T::residualPolicy_, T::jacobianPolicy_, auxStepperObj_)
  {}
};

// aux_stepper_t == void
template <typename T, typename stepper_t>
struct ImplicitStepperMixin<T, void, stepper_t> : T
{
  stepper_t stepperObj_;

  ImplicitStepperMixin() = delete;
  ImplicitStepperMixin(const ImplicitStepperMixin &) = default;
  ImplicitStepperMixin & operator=(const ImplicitStepperMixin &) = delete;
  ImplicitStepperMixin(ImplicitStepperMixin &&) = default;
  ImplicitStepperMixin & operator=(ImplicitStepperMixin &&) = delete;
  ~ImplicitStepperMixin() = default;

  template<typename T1, typename...Args>
  ImplicitStepperMixin(const T1 & romStateIn, Args && ...args)
    : T(romStateIn, std::forward<Args>(args)...),
      stepperObj_(romStateIn, T::fomCRef(),
		  T::residualPolicy_, T::jacobianPolicy_)
  {}
};

//---------------------------------------------------
// explicit stepper mixin
//---------------------------------------------------
template <typename T, typename stepper_t>
struct ExplicitStepperMixin : T
{
  stepper_t stepperObj_;

  ExplicitStepperMixin() = delete;
  ExplicitStepperMixin(const ExplicitStepperMixin &) = default;
  ExplicitStepperMixin & operator=(const ExplicitStepperMixin &) = delete;
  ExplicitStepperMixin(ExplicitStepperMixin &&) = default;
  ExplicitStepperMixin & operator=(ExplicitStepperMixin &&) = delete;
  ~ExplicitStepperMixin() = default;

  template<typename T1, typename...Args>
  ExplicitStepperMixin(const T1 & romStateIn,
		       Args && ...args)
    : T(romStateIn, std::forward<Args>(args)...),
      stepperObj_(romStateIn, T::fomCRef(), T::rhsPolicy_)
  {}
};

}}}
#endif  // ROM_IMPL_ROM_PROBLEM_MEMBERS_MIXINS_HPP_
