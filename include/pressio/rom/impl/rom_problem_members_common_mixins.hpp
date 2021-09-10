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
struct FomObjHolder;

template<class fom_system_t>
struct FomObjHolder<fom_system_t, false>
{
  std::reference_wrapper<const fom_system_t> fomObj_;

  FomObjHolder() = delete;
  FomObjHolder(const FomObjHolder &) = default;
  FomObjHolder & operator=(const FomObjHolder &) = delete;
  FomObjHolder(FomObjHolder &&) = default;
  FomObjHolder & operator=(FomObjHolder &&) = delete;
  ~FomObjHolder() = default;

  explicit FomObjHolder(const fom_system_t & fomObjIn) : fomObj_(fomObjIn){}

  const fom_system_t & fomCRef() const{ return fomObj_.get(); }
};

template <class T, class StepperType>
struct AddExplicitStepper : T
{
  StepperType stepperObj_;

  AddExplicitStepper() = delete;
  AddExplicitStepper(const AddExplicitStepper &) = default;
  AddExplicitStepper & operator=(const AddExplicitStepper &) = delete;
  AddExplicitStepper(AddExplicitStepper &&) = default;
  AddExplicitStepper & operator=(AddExplicitStepper &&) = delete;
  ~AddExplicitStepper() = default;

  template<class T1, class...Args>
  AddExplicitStepper(::pressio::ode::StepScheme name,
		       const T1 & romState,
		       Args && ...args)
    : T(name, romState, std::forward<Args>(args)...),
      stepperObj_(::pressio::ode::create_explicit_stepper(name, romState, T::romCRef()))
  {}
};

//---------------------------------------------------
template <class T, class StepperType>
struct AddImplicitStepper : T
{
  StepperType stepperObj_;

  AddImplicitStepper() = delete;
  AddImplicitStepper(const AddImplicitStepper &) = default;
  AddImplicitStepper & operator=(const AddImplicitStepper &) = delete;
  AddImplicitStepper(AddImplicitStepper &&) = default;
  AddImplicitStepper & operator=(AddImplicitStepper &&) = delete;
  ~AddImplicitStepper() = default;

  template<class T1, class...Args>
  AddImplicitStepper(::pressio::ode::StepScheme name,
		     const T1 & romState,
		     Args && ...args)
    : T(name, romState, std::forward<Args>(args)...),
      stepperObj_(::pressio::ode::create_implicit_stepper(name,
							  romState,
							  T::residualPolicy_,
							  T::jacobianPolicy_))
  {}
};

//---------------------------------------------------
template <class T, class StepperType>
struct AddImplicitArbStepper : T
{
  StepperType stepperObj_;

  AddImplicitArbStepper() = delete;
  AddImplicitArbStepper(const AddImplicitArbStepper &) = default;
  AddImplicitArbStepper & operator=(const AddImplicitArbStepper &) = delete;
  AddImplicitArbStepper(AddImplicitArbStepper &&) = default;
  AddImplicitArbStepper & operator=(AddImplicitArbStepper &&) = delete;
  ~AddImplicitArbStepper() = default;

  template<class T1, class...Args>
  AddImplicitArbStepper(::pressio::ode::StepScheme name,
			const T1 & romState,
			Args && ...args)
    : T(name, romState, std::forward<Args>(args)...),
      stepperObj_(romState, T::romCRef())
  {}
};

}}}
#endif  // ROM_IMPL_ROM_PROBLEM_MEMBERS_MIXINS_HPP_




// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template<class fom_system_t>
// struct FomObjHolder<fom_system_t, true>
// {
//   // when dealing with bindings for pressio4py, the fom_system_t
//   // is a C++ class in pressio4py that wraps the actual FOM python object.
//   // to construct this ROM problem, the Python code passes the
//   // python FOM object NOT a C++ object instantiated from fom_system_t.
//   // Therefore, ONLY when we deal with pressio4py, we create the fom obj
//   // instead of referencing it.
//   fom_system_t fomObj_;

//   FomObjHolder() = delete;
//   FomObjHolder(const FomObjHolder &) = default;
//   FomObjHolder & operator=(const FomObjHolder &) = delete;
//   FomObjHolder(FomObjHolder &&) = default;
//   FomObjHolder & operator=(FomObjHolder &&) = delete;
//   ~FomObjHolder() = default;

//   explicit FomObjHolder(pybind11::object pyFomObj): fomObj_(pyFomObj){}
//   const fom_system_t & fomCRef() const{ return fomObj_; }
// };
// #endif