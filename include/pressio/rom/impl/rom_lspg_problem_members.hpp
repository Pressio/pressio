/*
//@HEADER
// ************************************************************************
//
// rom_problem_members.hpp
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

#ifndef ROM_LSPG_IMPL_ROM_PROBLEM_MEMBERS_STEADY_HPP_
#define ROM_LSPG_IMPL_ROM_PROBLEM_MEMBERS_STEADY_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template <
  class T,
  class FomStateType,
  class FomStateReconstructorType,
  class ManagerStencilFomStatesType
  >
struct AddFomStatesManagerSteady : T
{
  const FomStateType  fomNominalState_;
  const FomStateReconstructorType fomStateReconstructor_;
  ManagerStencilFomStatesType fomStatesMngr_;

  AddFomStatesManagerSteady() = delete;
  AddFomStatesManagerSteady(const AddFomStatesManagerSteady &) = default;
  AddFomStatesManagerSteady & operator=(const AddFomStatesManagerSteady &) = delete;
  AddFomStatesManagerSteady(AddFomStatesManagerSteady &&) = default;
  AddFomStatesManagerSteady & operator=(AddFomStatesManagerSteady &&) = delete;
  ~AddFomStatesManagerSteady() = default;

  template<class T1, class T2, class T3, class T4>
  AddFomStatesManagerSteady(const T1 & fomObj,
			    const T2 & decoder,
			    const T3 & romStateIn,
			    const T4 & fomNominalStateNative)
    : T(fomObj),
      fomNominalState_(::pressio::ops::clone(fomNominalStateNative)),
      fomStateReconstructor_(fomNominalState_, decoder),
      fomStatesMngr_(fomStateReconstructor_, fomNominalState_)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(romStateIn);
  }
};

// ---------------------------------------------------------------
/* fom states manager
   need to store as many fom states as required by the scheme stencil
   - bdf1 needs 2
   - bdf2 needs 3
   - cn needs 2
*/

template<
  class ManagerStencilFomStatesType,
  class FomStateReconstructorType,
  class FomStateType>
ManagerStencilFomStatesType
create_manager_stencil_fom_states(::pressio::ode::SteppersE name,
				  const FomStateReconstructorType & fomRec,
				  const FomStateType & fomNomState)
{
  const auto tmp_b = ::pressio::ode::is_explicit_scheme(name);
  if (tmp_b){
    throw std::runtime_error("Unsteady LSPG prob members: enum must be implicit");
  }
  else{
    if (name == ::pressio::ode::SteppersE::BDF1){
      return ManagerStencilFomStatesType(fomRec,
					 {::pressio::ops::clone(fomNomState),
					  ::pressio::ops::clone(fomNomState)});
    }
    else if (name == ::pressio::ode::SteppersE::BDF2)
    {
      return ManagerStencilFomStatesType(fomRec,
					 {::pressio::ops::clone(fomNomState),
					  ::pressio::ops::clone(fomNomState),
					  ::pressio::ops::clone(fomNomState)});
    }
    else if (name == ::pressio::ode::SteppersE::CrankNicolson){
      return ManagerStencilFomStatesType(fomRec,
					 {::pressio::ops::clone(fomNomState),
					  ::pressio::ops::clone(fomNomState)}
					 );
    }
    else{
      throw std::runtime_error("Unsteady LSPG prob members: Invalid enum value");
    }
  }
}

template <
  class T,
  bool is_static,
  class FomStateType,
  class FomStateReconstructorType,
  class ManagerStencilFomStatesType
  >
struct AddFomStatesManagerUnsteady : T
{
  const FomStateType  fomNominalState_;
  const FomStateReconstructorType fomStateReconstructor_;
  ManagerStencilFomStatesType fomStatesMngr_;

  AddFomStatesManagerUnsteady() = delete;
  AddFomStatesManagerUnsteady(const AddFomStatesManagerUnsteady &) = default;
  AddFomStatesManagerUnsteady & operator=(const AddFomStatesManagerUnsteady &) = delete;
  AddFomStatesManagerUnsteady(AddFomStatesManagerUnsteady &&) = default;
  AddFomStatesManagerUnsteady & operator=(AddFomStatesManagerUnsteady &&) = delete;
  ~AddFomStatesManagerUnsteady() = default;

  template<
    class T1, class T2, class T3, class T4,
    bool _is_static = is_static,
    mpl::enable_if_t<!_is_static, int> = 0
    >
  AddFomStatesManagerUnsteady(::pressio::ode::SteppersE name,
			      const T1 & fomObj,
			      const T2 & decoder,
			      const T3 & romStateIn,
			      const T4 & fomNominalStateNative)
    : T(fomObj),
      fomNominalState_(::pressio::ops::clone(fomNominalStateNative)),
      fomStateReconstructor_(fomNominalState_, decoder),
      fomStatesMngr_(create_manager_stencil_fom_states<
		     ManagerStencilFomStatesType>(name,
						  fomStateReconstructor_,
						  fomNominalState_)
		     )
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(romStateIn);
  }

  template<
    class T1, class T2, class T3, class T4,
    bool _is_static = is_static,
    mpl::enable_if_t<_is_static, int> = 0
    >
  AddFomStatesManagerUnsteady(const T1 & fomObj,
			      const T2 & decoder,
			      const T3 & romStateIn,
			      const T4 & fomNominalStateNative)
    : T(fomObj),
      fomNominalState_(::pressio::ops::clone(fomNominalStateNative)),
      fomStateReconstructor_(fomNominalState_, decoder),
      fomStatesMngr_(fomStateReconstructor_, fomNominalState_)
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(romStateIn);
  }
};


template <class T, class r_pol_t, class j_pol_t>
struct AddDefaultPolicies : T
{
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  AddDefaultPolicies() = delete;
  AddDefaultPolicies(const AddDefaultPolicies &) = default;
  AddDefaultPolicies & operator=(const AddDefaultPolicies &) = delete;
  AddDefaultPolicies(AddDefaultPolicies &&) = default;
  AddDefaultPolicies & operator=(AddDefaultPolicies &&) = delete;
  ~AddDefaultPolicies() = default;

  template<class T1, class T2, class T3, class T4>
  AddDefaultPolicies(const T1 & romStateIn,
		     const T2 & fomObj,
		     T3 & decoder,
		     const T4 & fomNominalState)
    : T(fomObj, decoder, romStateIn, fomNominalState),
      residualPolicy_(T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(T::fomCRef(), T::fomStatesMngr_, decoder)
  {}

  template<class T1, class T2, class T3, class T4>
  AddDefaultPolicies(::pressio::ode::SteppersE name,
		     const T1 & romStateIn,
		     const T2 & fomObj,
		     T3 & decoder,
		     const T4 & fomNominalState)
    : T(name, fomObj, decoder, romStateIn, fomNominalState),
      residualPolicy_(T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(T::fomCRef(), T::fomStatesMngr_, decoder)
  {}
};


template <class T, class r_pol_t, class j_pol_t>
struct AddHypRedPolicies : T
{
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  AddHypRedPolicies() = delete;
  AddHypRedPolicies(const AddHypRedPolicies &) = default;
  AddHypRedPolicies & operator=(const AddHypRedPolicies &) = delete;
  AddHypRedPolicies(AddHypRedPolicies &&) = default;
  AddHypRedPolicies & operator=(AddHypRedPolicies &&) = delete;
  ~AddHypRedPolicies() = default;

  template<class T1, class T2, class T3, class T4, class T5>
  AddHypRedPolicies(::pressio::ode::SteppersE name,
		    const T1 & romStateIn,
		    const T2 & fomObj,
		    T3 & decoder,
		    const T4 & fomNominalState,
		    const T5 & combiner)
    : T(name, fomObj, decoder, romStateIn, fomNominalState),
      residualPolicy_(T::fomCRef(), T::fomStatesMngr_, combiner),
      jacobianPolicy_(T::fomCRef(), T::fomStatesMngr_, decoder, combiner)
  {}
};


template <class T, class UserProvidedFunctor_t, class r_pol_t, class j_pol_t>
struct AddSinglyDecoratedPolicies : T
{
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  AddSinglyDecoratedPolicies() = delete;
  AddSinglyDecoratedPolicies(const AddSinglyDecoratedPolicies &) = default;
  AddSinglyDecoratedPolicies & operator=(const AddSinglyDecoratedPolicies &) = delete;
  AddSinglyDecoratedPolicies(AddSinglyDecoratedPolicies &&) = default;
  AddSinglyDecoratedPolicies & operator=(AddSinglyDecoratedPolicies &&) = delete;
  ~AddSinglyDecoratedPolicies() = default;

  template<class T1, class T2, class T3, class T4>
  AddSinglyDecoratedPolicies(const T1 & romStateIn,
			     const T2 & fomObj,
			     T3 & decoder,
			     const T4 & fomNominalState,
			     const UserProvidedFunctor_t & userProvidedFunctor)
    : T(fomObj, decoder, romStateIn, fomNominalState),
      residualPolicy_(userProvidedFunctor, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(userProvidedFunctor, T::fomCRef(), T::fomStatesMngr_, decoder)
  {}

  template<class T1, class T2, class T3, class T4>
  AddSinglyDecoratedPolicies(::pressio::ode::SteppersE name,
			     const T1 & romStateIn,
			     const T2 & fomObj,
			     T3 & decoder,
			     const T4 & fomNominalState,
			     const UserProvidedFunctor_t & userProvidedFunctor)
    : T(name, fomObj, decoder, romStateIn, fomNominalState),
      residualPolicy_(userProvidedFunctor, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(userProvidedFunctor, T::fomCRef(), T::fomStatesMngr_, decoder)
  {}
};

template <
  class T,
  class UserProvidedFunctor1_t,
  class UserProvidedFunctor2_t,
  class r_pol_t,
  class j_pol_t
  >
struct AddDoublyDecoratedPolicies : T
{
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  AddDoublyDecoratedPolicies() = delete;
  AddDoublyDecoratedPolicies(const AddDoublyDecoratedPolicies &) = default;
  AddDoublyDecoratedPolicies & operator=(const AddDoublyDecoratedPolicies &) = delete;
  AddDoublyDecoratedPolicies(AddDoublyDecoratedPolicies &&) = default;
  AddDoublyDecoratedPolicies & operator=(AddDoublyDecoratedPolicies &&) = delete;
  ~AddDoublyDecoratedPolicies() = default;

  template<class T1, class T2, class T3, class T4>
  AddDoublyDecoratedPolicies(const T1 & romStateIn,
			     const T2 & fomObj,
			     T3 & decoder,
			     const T4 & fomNominalState,
			     const UserProvidedFunctor1_t & functor1,
			     const UserProvidedFunctor2_t & functor2)
    : T(fomObj, decoder, romStateIn, fomNominalState),
      residualPolicy_(functor1, functor2, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(functor1, functor2, T::fomCRef(), T::fomStatesMngr_, decoder)
  {}

  template<class T1, class T2, class T3, class T4>
  AddDoublyDecoratedPolicies(::pressio::ode::SteppersE name,
			     const T1 & romStateIn,
			     const T2 & fomObj,
			     T3 & decoder,
			     const T4 & fomNominalState,
			     const UserProvidedFunctor1_t & functor1,
			     const UserProvidedFunctor2_t & functor2)
    : T(name, fomObj, decoder, romStateIn, fomNominalState),
      residualPolicy_(functor1, functor2, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(functor1, functor2, T::fomCRef(), T::fomStatesMngr_, decoder)
  {}
};

template <class T, class rom_sys_t>
struct AddDefaultDiscreteTimeSystem : T
{
  const rom_sys_t romSys_;

  AddDefaultDiscreteTimeSystem() = delete;
  AddDefaultDiscreteTimeSystem(const AddDefaultDiscreteTimeSystem &) = default;
  AddDefaultDiscreteTimeSystem & operator=(const AddDefaultDiscreteTimeSystem &) = delete;
  AddDefaultDiscreteTimeSystem(AddDefaultDiscreteTimeSystem &&) = default;
  AddDefaultDiscreteTimeSystem & operator=(AddDefaultDiscreteTimeSystem &&) = delete;
  ~AddDefaultDiscreteTimeSystem() = default;

  template<class T1, class T2, class T3, class T4>
  AddDefaultDiscreteTimeSystem(::pressio::ode::SteppersE name,
			       const T1 & romStateIn,
			       const T2 & fomObj,
			       T3 & decoder,
			       const T4 & fomNominalState)
    : T(fomObj, decoder, romStateIn, fomNominalState),
      romSys_(T::fomCRef(), decoder, T::fomStatesMngr_)
  {}

  const rom_sys_t & romCRef() const{ return romSys_; }
};

template <class T, class UserProvidedFunctor_t, class rom_sys_t>
struct AddSinglyDecoratedDiscreteTimeSystem : T
{
  const rom_sys_t romSys_;

  AddSinglyDecoratedDiscreteTimeSystem() = delete;
  AddSinglyDecoratedDiscreteTimeSystem(const AddSinglyDecoratedDiscreteTimeSystem &) = default;
  AddSinglyDecoratedDiscreteTimeSystem & operator=(const AddSinglyDecoratedDiscreteTimeSystem &) = delete;
  AddSinglyDecoratedDiscreteTimeSystem(AddSinglyDecoratedDiscreteTimeSystem &&) = default;
  AddSinglyDecoratedDiscreteTimeSystem & operator=(AddSinglyDecoratedDiscreteTimeSystem &&) = delete;
  ~AddSinglyDecoratedDiscreteTimeSystem() = default;

  template<class T1, class T2, class T3, class T4>
  AddSinglyDecoratedDiscreteTimeSystem(::pressio::ode::SteppersE name,
				       const T1 & romStateIn,
				       const T2 & fomObj,
				       T3 & decoder,
				       const T4 & fomNominalState,
				       const UserProvidedFunctor_t & userProvidedFunctor)
    : T(fomObj, decoder, romStateIn, fomNominalState),
      romSys_(userProvidedFunctor, T::fomCRef(), decoder, T::fomStatesMngr_)
  {}

  const rom_sys_t & romCRef() const{ return romSys_; }
};


template <
  class T,
  class UserProvidedFunctor1_t,
  class UserProvidedFunctor2_t,
  class rom_sys_t>
struct AddDoublyDecoratedDiscreteTimeSystem : T
{
  const rom_sys_t romSys_;

  AddDoublyDecoratedDiscreteTimeSystem() = delete;
  AddDoublyDecoratedDiscreteTimeSystem(const AddDoublyDecoratedDiscreteTimeSystem &) = default;
  AddDoublyDecoratedDiscreteTimeSystem & operator=(const AddDoublyDecoratedDiscreteTimeSystem &) = delete;
  AddDoublyDecoratedDiscreteTimeSystem(AddDoublyDecoratedDiscreteTimeSystem &&) = default;
  AddDoublyDecoratedDiscreteTimeSystem & operator=(AddDoublyDecoratedDiscreteTimeSystem &&) = delete;
  ~AddDoublyDecoratedDiscreteTimeSystem() = default;

  template<class T1, class T2, class T3, class T4>
  AddDoublyDecoratedDiscreteTimeSystem(::pressio::ode::SteppersE name,
				       const T1 & romStateIn,
				       const T2 & fomObj,
				       T3 & decoder,
				       const T4 & fomNominalState,
				       const UserProvidedFunctor1_t & userProvidedFunctor1,
				       const UserProvidedFunctor2_t & userProvidedFunctor2)
    : T(fomObj, decoder, romStateIn, fomNominalState),
      romSys_(userProvidedFunctor1, userProvidedFunctor2,
	      T::fomCRef(), decoder, T::fomStatesMngr_)
  {}

  const rom_sys_t & romCRef() const{ return romSys_; }
};


}}}}
#endif  // ROM_LSPG_IMPL_ROM_PROBLEM_MEMBERS_HPP_
