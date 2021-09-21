/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_problem_members.hpp
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

#ifndef ROM_IMPL_ROM_GALERKIN_PROBLEM_MEMBERS_HPP_
#define ROM_IMPL_ROM_GALERKIN_PROBLEM_MEMBERS_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{


// ---------------------------------------------------------------
/* fom states manager

   - for Galerkin the FOM states are only needed when we query the FOM velocity.

   - For explicit time stepping, we only need to store one fom state that
   we reconstruct every time we need to compute the FOM velocity

   - For implicit time stepping, we might need to store more than one
   depending on how many FOM velocity evaluations are needed.
   BDF1 and BDF2 need FOM velocity at n+1, so one FOM state.
   CrankNicolson needs two evaluations of the FOM velocity
   at n and n+1, so we need two FOM states.
*/

template<
  class ManagerStencilFomStatesType,
  class FomStateReconstructorType,
  class FomStateType>
ManagerStencilFomStatesType
create_manager_stencil_fom_states(::pressio::ode::StepScheme name,
				  const FomStateReconstructorType & fomRec,
				  const FomStateType & fomNomState)
{
  const auto tmp_b = ::pressio::ode::is_explicit_scheme(name);
  if (tmp_b){
    // use initializer_list here!
    return ManagerStencilFomStatesType(fomRec, {::pressio::ops::clone(fomNomState)});
  }
  else{
    if (name == ::pressio::ode::StepScheme::BDF1 or name == ::pressio::ode::StepScheme::BDF2){
      return ManagerStencilFomStatesType(fomRec, {::pressio::ops::clone(fomNomState)});
    }
    else if (name == ::pressio::ode::StepScheme::CrankNicolson){
      return ManagerStencilFomStatesType(fomRec,
					 {::pressio::ops::clone(fomNomState),
					  ::pressio::ops::clone(fomNomState)}
					 );
    }
    else{
      throw std::runtime_error("Galerkin prob members: Invalid enum value");
    }
  }
}

template <
  class T,
  bool is_cont_time,
  class FomStateType,
  class FomStateReconstructorType,
  class ManagerStencilFomStatesType
  >
struct AddFomStatesManager : T
{
  const FomStateType  fomNominalState_;
  const FomStateReconstructorType fomStateReconstructor_;
  ManagerStencilFomStatesType fomStatesMngr_;

  AddFomStatesManager() = delete;
  AddFomStatesManager(const AddFomStatesManager &) = default;
  AddFomStatesManager & operator=(const AddFomStatesManager &) = delete;
  AddFomStatesManager(AddFomStatesManager &&) = default;
  AddFomStatesManager & operator=(AddFomStatesManager &&) = delete;
  ~AddFomStatesManager() = default;

  template<
    class T1, class T2, class T3, class T4,
    bool _is_cont_time = is_cont_time,
    mpl::enable_if_t<_is_cont_time, int> = 0
    >
  AddFomStatesManager(::pressio::ode::StepScheme name,
		      const T1 & fomObj,
		      const T2 & decoder,
		      const T3 & romStateIn,
		      const T4 & fomNominalStateNative)
    : T(fomObj),
      fomNominalState_(::pressio::ops::clone(fomNominalStateNative)),
      fomStateReconstructor_(fomNominalState_, decoder),
      fomStatesMngr_(create_manager_stencil_fom_states<
		     ManagerStencilFomStatesType>(name, fomStateReconstructor_, fomNominalState_)
		     )
  {
    // reconstruct current fom state so that we have something
    // consisten with the current romState
    fomStatesMngr_.reconstructCurrentFomState(romStateIn);
  }

  template<
    class T1, class T2, class T3, class T4,
    bool _is_cont_time = is_cont_time,
    mpl::enable_if_t<!_is_cont_time, int> = 0
    >
  AddFomStatesManager(::pressio::ode::StepScheme name,
		      const T1 & fomObj,
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

//---------------------------------------------------
template <class T, class projector_t>
struct AddProjector : T
{
  const projector_t projector_;

  AddProjector() = delete;
  AddProjector(const AddProjector &) = default;
  AddProjector & operator=(const AddProjector &) = delete;
  AddProjector(AddProjector &&) = default;
  AddProjector & operator=(AddProjector &&) = delete;
  ~AddProjector() = default;

  template<class T1, class T2, class T3, class T4>
  AddProjector(::pressio::ode::StepScheme name,
		 const T1 & romStateIn,
		 const T2 & fomObj,
		 const T3 & decoder,
		 const T4 & fomNominalState)
    : T(name, fomObj, decoder, romStateIn, fomNominalState),
      projector_(decoder)
  {}

  template<class T1, class T2, class T3, class T4, class T5, class ...Args>
  AddProjector(::pressio::ode::StepScheme name,
	       const T1 & romStateIn,
	       const T2 & fomObj,
	       const T3 & decoder,
	       const T4 & fomNominalState,
	       const T5 & projector,
	       Args && ... args)
    : T(name, fomObj, decoder, romStateIn, fomNominalState, std::forward<Args>(args)...),
      projector_(projector)
  {}
};

template <class T, class masker_t>
struct AddMasker : T
{
  const masker_t masker_;

  AddMasker() = delete;
  AddMasker(const AddMasker &) = default;
  AddMasker & operator=(const AddMasker &) = delete;
  AddMasker(AddMasker &&) = default;
  AddMasker & operator=(AddMasker &&) = delete;
  ~AddMasker() = default;

  template<class T1, class T2, class T3, class T4, class T5, class ...Args>
  AddMasker(::pressio::ode::StepScheme name,
	      const T1 & fomObj,
	      const T2 & decoder,
	      const T3 & romStateIn,
	      const T4 & fomNominalState,
	      const T5 & masker,
	      Args && ... args)
    : T(name, fomObj, decoder, romStateIn, fomNominalState, std::forward<Args>(args)...),
      masker_(masker)
  {}
};

template <class T, class rom_sys_t>
struct AddDefaultExplicitSystem : T
{
  const rom_sys_t romSys_;

  AddDefaultExplicitSystem() = delete;
  AddDefaultExplicitSystem(const AddDefaultExplicitSystem &) = default;
  AddDefaultExplicitSystem & operator=(const AddDefaultExplicitSystem &) = delete;
  AddDefaultExplicitSystem(AddDefaultExplicitSystem &&) = default;
  AddDefaultExplicitSystem & operator=(AddDefaultExplicitSystem &&) = delete;
  ~AddDefaultExplicitSystem() = default;

  template<class T1, typename ...Args>
  AddDefaultExplicitSystem(::pressio::ode::StepScheme name,
			   const T1 & romStateIn,
			   Args && ...args)
    : T(name, romStateIn, std::forward<Args>(args)...),
      romSys_(romStateIn, T::projector_, T::fomCRef(), T::fomStatesMngr_)
  {}

  const rom_sys_t & romCRef() const{ return romSys_; }
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
  AddDefaultDiscreteTimeSystem(::pressio::ode::StepScheme name,
			       const T1 & romStateIn,
			       const T2 & fomObj,
			       const T3 & decoder,
			       const T4 & fomNominalState)
    : T(name, romStateIn, fomObj, decoder, fomNominalState),
      romSys_(romStateIn, T::projector_, T::fomCRef(), decoder, T::fomStatesMngr_)
  {}

  const rom_sys_t & romCRef() const{ return romSys_; }
};

template <class T, class rom_sys_t>
struct AddHypRedDiscreteTimeSystem : T
{
  const rom_sys_t romSys_;

  AddHypRedDiscreteTimeSystem() = delete;
  AddHypRedDiscreteTimeSystem(const AddHypRedDiscreteTimeSystem &) = default;
  AddHypRedDiscreteTimeSystem & operator=(const AddHypRedDiscreteTimeSystem &) = delete;
  AddHypRedDiscreteTimeSystem(AddHypRedDiscreteTimeSystem &&) = default;
  AddHypRedDiscreteTimeSystem & operator=(AddHypRedDiscreteTimeSystem &&) = delete;
  ~AddHypRedDiscreteTimeSystem() = default;

  template<class T1, class T2, class T3, class T4, class T5>
  AddHypRedDiscreteTimeSystem(::pressio::ode::StepScheme name,
			      const T1 & romStateIn,
			      const T2 & fomObj,
			      const T3 & decoder,
			      const T4 & fomNominalState,
			      const T5 & projector)
    : T(name, romStateIn, fomObj, decoder, fomNominalState, projector),
      romSys_(romStateIn, T::projector_, T::fomCRef(), decoder, T::fomStatesMngr_)
  {}

  const rom_sys_t & romCRef() const{ return romSys_; }
};

template <class T, class rom_sys_t>
struct AddMaskedDiscreteTimeSystem : T
{
  const rom_sys_t romSys_;

  AddMaskedDiscreteTimeSystem() = delete;
  AddMaskedDiscreteTimeSystem(const AddMaskedDiscreteTimeSystem &) = default;
  AddMaskedDiscreteTimeSystem & operator=(const AddMaskedDiscreteTimeSystem &) = delete;
  AddMaskedDiscreteTimeSystem(AddMaskedDiscreteTimeSystem &&) = default;
  AddMaskedDiscreteTimeSystem & operator=(AddMaskedDiscreteTimeSystem &&) = delete;
  ~AddMaskedDiscreteTimeSystem() = default;

  template<class T1, class T2, class T3, class T4, class T5, class T6>
  AddMaskedDiscreteTimeSystem(::pressio::ode::StepScheme name,
			      const T1 & romStateIn,
			      const T2 & fomObj,
			      const T3 & decoder,
			      const T4 & fomNominalState,
			      const T5 & projector,
			      const T6 & masker)
    : T(name, romStateIn, fomObj, decoder, fomNominalState, projector, masker),
      romSys_(romStateIn, T::projector_, T::masker_, T::fomCRef(), decoder, T::fomStatesMngr_)
  {}

  const rom_sys_t & romCRef() const{ return romSys_; }
};

template <class T, class rom_sys_t>
struct AddMaskedVeloExplicitSystem : T
{
  const rom_sys_t romSys_;

  AddMaskedVeloExplicitSystem() = delete;
  AddMaskedVeloExplicitSystem(const AddMaskedVeloExplicitSystem &) = default;
  AddMaskedVeloExplicitSystem & operator=(const AddMaskedVeloExplicitSystem &) = delete;
  AddMaskedVeloExplicitSystem(AddMaskedVeloExplicitSystem &&) = default;
  AddMaskedVeloExplicitSystem & operator=(AddMaskedVeloExplicitSystem &&) = delete;
  ~AddMaskedVeloExplicitSystem() = default;

  template<class T1, typename ...Args>
  AddMaskedVeloExplicitSystem(::pressio::ode::StepScheme name,
			      const T1 & romStateIn,
			      Args && ...args)
    : T(name, romStateIn, std::forward<Args>(args)...),
      romSys_(romStateIn, T::projector_, T::masker_, T::fomCRef(), T::fomStatesMngr_)
  {}

  const rom_sys_t & romCRef() const{ return romSys_; }
};

template <class T, class rom_sys_t>
using AddHypRedVeloExplicitSystem = AddDefaultExplicitSystem<T, rom_sys_t>;

//---------------------------------------------------
// implicit policies
//---------------------------------------------------
template <
  class T,
  bool is_default, bool is_hypredvelo, bool is_masked,
  class ... Args>
struct ImplicitPoliciesMixin;

// specialize for default
template <class T, class r_pol_t, class j_pol_t>
struct ImplicitPoliciesMixin<T, true, false, false, r_pol_t, j_pol_t> : T
{
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  ImplicitPoliciesMixin() = delete;
  ImplicitPoliciesMixin(const ImplicitPoliciesMixin &) = default;
  ImplicitPoliciesMixin & operator=(const ImplicitPoliciesMixin &) = delete;
  ImplicitPoliciesMixin(ImplicitPoliciesMixin &&) = default;
  ImplicitPoliciesMixin & operator=(ImplicitPoliciesMixin &&) = delete;
  ~ImplicitPoliciesMixin() = default;

  template<class T1, class T2, class T3, class T4, class ...Args>
  ImplicitPoliciesMixin(::pressio::ode::StepScheme name,
			const T1 & romStateIn,
			const T2 & fomObj,
			const T3 & decoder,
			const T4 & fomNominalState,
			Args && ... args)
    : T(name, romStateIn, fomObj, decoder, fomNominalState, std::forward<Args>(args)...),
      residualPolicy_(::pressio::ops::extent(romStateIn,0),
		      T::projector_, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(::pressio::ops::extent(romStateIn,0),
		      ::pressio::ops::extent(romStateIn,0),
		      T::projector_, T::fomCRef(), T::fomStatesMngr_, decoder)
  {}
};

// specialize for hyp-redu velo
template <class T, class r_pol_t, class j_pol_t>
struct ImplicitPoliciesMixin<T, false, true, false, r_pol_t, j_pol_t> : T
{
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  ImplicitPoliciesMixin() = delete;
  ImplicitPoliciesMixin(const ImplicitPoliciesMixin &) = default;
  ImplicitPoliciesMixin & operator=(const ImplicitPoliciesMixin &) = delete;
  ImplicitPoliciesMixin(ImplicitPoliciesMixin &&) = default;
  ImplicitPoliciesMixin & operator=(ImplicitPoliciesMixin &&) = delete;
  ~ImplicitPoliciesMixin() = default;

  template<class T1, class T2, class T3, class T4, class ...Args>
  ImplicitPoliciesMixin(::pressio::ode::StepScheme name,
			const T1 & romStateIn,
			const T2 & fomObj,
			const T3 & decoder,
			const T4 & fomNominalState,
			Args && ... args)
    : T(name, romStateIn, fomObj, decoder, fomNominalState, std::forward<Args>(args)...),
      residualPolicy_(::pressio::ops::extent(romStateIn,0),
		      T::projector_, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(::pressio::ops::extent(romStateIn,0),
		      ::pressio::ops::extent(romStateIn,0),
		      T::projector_, T::fomCRef(), T::fomStatesMngr_, decoder)
  {}
};

// specialize for masked velo
template <class T, class r_pol_t, class j_pol_t>
struct ImplicitPoliciesMixin<T, false, false, true, r_pol_t, j_pol_t> : T
{
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  ImplicitPoliciesMixin() = delete;
  ImplicitPoliciesMixin(const ImplicitPoliciesMixin &) = default;
  ImplicitPoliciesMixin & operator=(const ImplicitPoliciesMixin &) = delete;
  ImplicitPoliciesMixin(ImplicitPoliciesMixin &&) = default;
  ImplicitPoliciesMixin & operator=(ImplicitPoliciesMixin &&) = delete;
  ~ImplicitPoliciesMixin() = default;

  template<class T1, class T2, class T3, class T4, class ...Args>
  ImplicitPoliciesMixin(::pressio::ode::StepScheme name,
			const T1 & romStateIn,
			const T2 & fomObj,
			const T3 & decoder,
			const T4 & fomNominalState,
			Args && ... args)
    : T(name, romStateIn, fomObj, decoder, fomNominalState, std::forward<Args>(args)...),
      residualPolicy_(::pressio::ops::extent(romStateIn,0),
		      T::projector_, T::masker_, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(::pressio::ops::extent(romStateIn,0),
		      ::pressio::ops::extent(romStateIn,0),
		      T::projector_, T::masker_, T::fomCRef(), T::fomStatesMngr_, decoder)
  {}
};

// aliases to make things easier
template <class T, typename ...Args>
using AddDefaultImplicitPolicies = ImplicitPoliciesMixin<T, true, false, false, Args...>;

template <class T, typename ...Args>
using AddHypRedVeloImplicitPolicies = ImplicitPoliciesMixin<T, false, true, false, Args...>;

template <class T, typename ...Args>
using AddMaskedVeloImplicitPolicies = ImplicitPoliciesMixin<T, false, false, true, Args...>;

}}}}
#endif  // ROM_IMPL_ROM_GALERKIN_PROBLEM_MEMBERS_HPP_
