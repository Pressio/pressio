/*
//@HEADER
// ************************************************************************
//
// rom_lspg_problem_members.hpp
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

#ifndef ROM_LSPG_IMPL_ROM_LSPG_PROBLEM_MEMBERS_HPP_
#define ROM_LSPG_IMPL_ROM_LSPG_PROBLEM_MEMBERS_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

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

  FomObjMixin(const fom_system_t & fomObjIn) : fomObj_(fomObjIn){}

  const fom_system_t & fomCRef() const{ return fomObj_.get(); }
};

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

  FomObjMixin(fom_system_t fomObjIn) : fomObj_(fomObjIn){}

  const fom_system_t & fomCRef() const{ return fomObj_; }
};

//---------------------------------------------------
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
//---------------------------------------------------
template <
  typename T,
  bool def, bool masked, bool prec, bool hypred,
  typename ... Args
  >
struct PoliciesMixin;

// specialize for default
template <typename T, typename ops_t, typename r_pol_t, typename j_pol_t>
struct PoliciesMixin<
  T, true, false, false, false, ops_t, r_pol_t, j_pol_t
  > : T
{
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  PoliciesMixin() = delete;
  PoliciesMixin(const PoliciesMixin &) = default;
  PoliciesMixin & operator=(const PoliciesMixin &) = delete;
  PoliciesMixin(PoliciesMixin &&) = default;
  PoliciesMixin & operator=(PoliciesMixin &&) = delete;
  ~PoliciesMixin() = default;

  template<
    typename T1, typename T2, typename T3, typename T4, typename _ops_t = ops_t,
    mpl::enable_if_t<std::is_void<_ops_t>::value, int > = 0
    >
  PoliciesMixin(const T1 & romStateIn,
		const T2 & fomObj,
		const T3 & decoder,
		const T4 & fomNominalStateNative)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative),
      residualPolicy_(T::fomStatesMngr_),
      jacobianPolicy_(T::fomStatesMngr_, decoder)
  {}

  template<
    typename T1, typename T2, typename T3, typename T4, typename _ops_t = ops_t,
    mpl::enable_if_t<!std::is_void<_ops_t>::value, int > = 0
    >
  PoliciesMixin(const T1 & romStateIn,
		const T2 & fomObj,
		const T3 & decoder,
		const T4 & fomNominalStateNative,
		const _ops_t & ops)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative, ops),
      residualPolicy_(T::fomStatesMngr_, ops),
      jacobianPolicy_(T::fomStatesMngr_, decoder, ops)
  {}
};

// specialize for masked
template <typename T, typename ops_t, typename r_pol_t, typename j_pol_t>
struct PoliciesMixin<
  T, false, true, false, false, ops_t, r_pol_t, j_pol_t
  > : T
{
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  PoliciesMixin() = delete;
  PoliciesMixin(const PoliciesMixin &) = default;
  PoliciesMixin & operator=(const PoliciesMixin &) = delete;
  PoliciesMixin(PoliciesMixin &&) = default;
  PoliciesMixin & operator=(PoliciesMixin &&) = delete;
  ~PoliciesMixin() = default;

  template<
    typename T1, typename T2, typename T3, typename T4, typename T5,
    typename _ops_t = ops_t,
    mpl::enable_if_t<std::is_void<_ops_t>::value, int > = 0
    >
  PoliciesMixin(const T1 & romStateIn,
		const T2 & fomObj,
		const T3 & decoder,
		const T4 & fomNominalStateNative,
		const T5 & masker)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative),
      residualPolicy_(masker, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(masker, T::fomCRef(), T::fomStatesMngr_, decoder)
  {}
};

// specialize for precond
template <typename T, typename ops_t, typename r_pol_t, typename j_pol_t>
struct PoliciesMixin<
  T, false, false, true, false, ops_t, r_pol_t, j_pol_t
  > : T
{
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  PoliciesMixin() = delete;
  PoliciesMixin(const PoliciesMixin &) = default;
  PoliciesMixin & operator=(const PoliciesMixin &) = delete;
  PoliciesMixin(PoliciesMixin &&) = default;
  PoliciesMixin & operator=(PoliciesMixin &&) = delete;
  ~PoliciesMixin() = default;

  template<
    typename T1, typename T2, typename T3, typename T4, typename T5,
    typename _ops_t = ops_t,
    mpl::enable_if_t<std::is_void<_ops_t>::value, int > = 0
    >
  PoliciesMixin(const T1 & romStateIn,
		const T2 & fomObj,
		const T3 & decoder,
		const T4 & fomNominalStateNative,
		const T5 & preconditioner)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative),
      residualPolicy_(preconditioner, T::fomStatesMngr_),
      jacobianPolicy_(preconditioner, T::fomStatesMngr_, decoder)
  {}
};

// specialize for hyp-red with void stencil-to-sample mapping
template <typename T, typename ops_t, typename r_pol_t, typename j_pol_t>
struct PoliciesMixin<
  T, false, false, false, true, ops_t, r_pol_t, j_pol_t, void
  > : T
{
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  PoliciesMixin() = delete;
  PoliciesMixin(const PoliciesMixin &) = default;
  PoliciesMixin & operator=(const PoliciesMixin &) = delete;
  PoliciesMixin(PoliciesMixin &&) = default;
  PoliciesMixin & operator=(PoliciesMixin &&) = delete;
  ~PoliciesMixin() = default;

  template<
    typename T1, typename T2, typename T3, typename T4,
    typename _ops_t = ops_t,
    mpl::enable_if_t<std::is_void<_ops_t>::value, int > = 0
    >
  PoliciesMixin(const T1 & romStateIn,
		const T2 & fomObj,
		const T3 & decoder,
		const T4 & fomNominalStateNative)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative),
      residualPolicy_(T::fomStatesMngr_),
      jacobianPolicy_(T::fomStatesMngr_, decoder)
  {}
};

// specialize for hyp-red with nonvoid stencil-to-sample mapping
template <typename T, typename ops_t, typename r_pol_t, typename j_pol_t, typename sTos_t>
struct PoliciesMixin<
  T, false, false, false, true, ops_t, r_pol_t, j_pol_t, sTos_t
  > : T
{
  sTos_t  meshToStencilMapper_;
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  PoliciesMixin() = delete;
  PoliciesMixin(const PoliciesMixin &) = default;
  PoliciesMixin & operator=(const PoliciesMixin &) = delete;
  PoliciesMixin(PoliciesMixin &&) = default;
  PoliciesMixin & operator=(PoliciesMixin &&) = delete;
  ~PoliciesMixin() = default;

  template<
    typename T1, typename T2, typename T3, typename T4, typename T5,
    typename _ops_t = ops_t,
    mpl::enable_if_t<std::is_void<_ops_t>::value, int > = 0
    >
  PoliciesMixin(const T1 & romStateIn,
		const T2 & fomObj,
		const T3 & decoder,
		const T4 & fomNominalStateNative,
		const T5 & meshToStencilMapper)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative),
      meshToStencilMapper_(meshToStencilMapper),
      residualPolicy_(T::fomStatesMngr_, meshToStencilMapper_),
      jacobianPolicy_(T::fomStatesMngr_, decoder, meshToStencilMapper_)
  {}
};

// aliases to make things easier
template <typename T, typename ...Args>
using DefaultPoliciesMixin = PoliciesMixin<T, true, false, false, false, Args...>;

template <typename T, typename ...Args>
using MaskedPoliciesMixin = PoliciesMixin<T, false, true, false, false, Args...>;

template <typename T, typename ...Args>
using PrecondPoliciesMixin = PoliciesMixin<T, false, false, true, false, Args...>;

template <typename T, typename ...Args>
using HypRedPoliciesMixin = PoliciesMixin<T, false, false, false, true, Args...>;



//---------------------------------------------------
// stepper mixin
//---------------------------------------------------
// aux_stepper_t is valid
template <typename T, typename aux_stepper_t, typename stepper_t>
struct StepperMixin : T
{
  aux_stepper_t auxStepperObj_;
  stepper_t stepperObj_;

  StepperMixin() = delete;
  StepperMixin(const StepperMixin &) = default;
  StepperMixin & operator=(const StepperMixin &) = delete;
  StepperMixin(StepperMixin &&) = default;
  StepperMixin & operator=(StepperMixin &&) = delete;
  ~StepperMixin() = default;

  template<typename T1, typename...Args>
  StepperMixin(const T1 & romStateIn,
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
struct StepperMixin<T, void, stepper_t> : T
{
  stepper_t stepperObj_;

  StepperMixin() = delete;
  StepperMixin(const StepperMixin &) = default;
  StepperMixin & operator=(const StepperMixin &) = delete;
  StepperMixin(StepperMixin &&) = default;
  StepperMixin & operator=(StepperMixin &&) = delete;
  ~StepperMixin() = default;

  template<typename T1, typename...Args>
  StepperMixin(const T1 & romStateIn,
	       Args && ...args)
    : T(romStateIn, std::forward<Args>(args)...),
      stepperObj_(romStateIn, T::fomCRef(),
		  T::residualPolicy_, T::jacobianPolicy_)
  {}
};


//---------------------------------------------------
// system mixin (used for steady LSPG)
//---------------------------------------------------
template <typename T, typename system_t>
struct SystemMixin : T
{
  system_t systemObj_;

  SystemMixin() = delete;
  SystemMixin(const SystemMixin &) = default;
  SystemMixin & operator=(const SystemMixin &) = delete;
  SystemMixin(SystemMixin &&) = default;
  SystemMixin & operator=(SystemMixin &&) = delete;
  ~SystemMixin() = default;

  template<typename...Args>
  SystemMixin(Args && ...args)
    : T(std::forward<Args>(args)...),
      systemObj_(T::fomCRef(), T::residualPolicy_, T::jacobianPolicy_)
  {}
};

}}}}
#endif  // ROM_LSPG_IMPL_ROM_LSPG_PROBLEM_MEMBERS_HPP_
