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

#ifndef ROM_GALERKIN_IMPL_ROM_PROBLEM_MEMBERS_HPP_
#define ROM_GALERKIN_IMPL_ROM_PROBLEM_MEMBERS_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <class T, class projector_t>
struct ProjectorMixin : T
{
  const projector_t projector_;

  ProjectorMixin() = delete;
  ProjectorMixin(const ProjectorMixin &) = default;
  ProjectorMixin & operator=(const ProjectorMixin &) = delete;
  ProjectorMixin(ProjectorMixin &&) = default;
  ProjectorMixin & operator=(ProjectorMixin &&) = delete;
  ~ProjectorMixin() = default;

  template<class T1, class T2, class T3, class T4>
  ProjectorMixin(const T1 & romStateIn,
		 const T2 & fomObj,
		 const T3 & decoder,
		 const T4 & fomNominalState)
    : T(fomObj, decoder, romStateIn, fomNominalState),
      projector_(decoder)
  {}

  template<class T1, class T2, class T3, class T4, class T5, class ...Args>
  ProjectorMixin(const T1 & romStateIn,
		 const T2 & fomObj,
		 const T3 & decoder,
		 const T4 & fomNominalState,
		 const T5 & projector,
		 Args && ... args)
    : T(fomObj, decoder, romStateIn, fomNominalState, std::forward<Args>(args)...),
      projector_(projector)
  {}
};

template <class T, class masker_t>
struct MaskerMixin : T
{
  const masker_t masker_;

  MaskerMixin() = delete;
  MaskerMixin(const MaskerMixin &) = default;
  MaskerMixin & operator=(const MaskerMixin &) = delete;
  MaskerMixin(MaskerMixin &&) = default;
  MaskerMixin & operator=(MaskerMixin &&) = delete;
  ~MaskerMixin() = default;

  template<class T1, class T2, class T3, class T4, class T5, class ...Args>
  MaskerMixin(const T1 & fomObj,
	      const T2 & decoder,
	      const T3 & romStateIn,
	      const T4 & fomNominalState,
	      const T5 & masker,
	      Args && ... args)
    : T(fomObj, decoder, romStateIn, fomNominalState, std::forward<Args>(args)...),
      masker_(masker)
  {}
};

template <class T, class rom_sys_t>
struct DefaultDiscreteTimeSystemMixin : T
{
  const rom_sys_t romSys_;

  DefaultDiscreteTimeSystemMixin() = delete;
  DefaultDiscreteTimeSystemMixin(const DefaultDiscreteTimeSystemMixin &) = default;
  DefaultDiscreteTimeSystemMixin & operator=(const DefaultDiscreteTimeSystemMixin &) = delete;
  DefaultDiscreteTimeSystemMixin(DefaultDiscreteTimeSystemMixin &&) = default;
  DefaultDiscreteTimeSystemMixin & operator=(DefaultDiscreteTimeSystemMixin &&) = delete;
  ~DefaultDiscreteTimeSystemMixin() = default;

  template<class T1, class T2, class T3, class T4>
  DefaultDiscreteTimeSystemMixin(const T1 & romStateIn,
				 const T2 & fomObj,
				 const T3 & decoder,
				 const T4 & fomNominalState)
    : T(romStateIn, fomObj, decoder, fomNominalState),
      romSys_(romStateIn, T::projector_, T::fomCRef(), decoder, T::fomStatesMngr_)
  {}

  const rom_sys_t & romCRef() const{ return romSys_; }
};

template <class T, class rom_sys_t>
struct HypRedDiscreteTimeSystemMixin : T
{
  const rom_sys_t romSys_;

  HypRedDiscreteTimeSystemMixin() = delete;
  HypRedDiscreteTimeSystemMixin(const HypRedDiscreteTimeSystemMixin &) = default;
  HypRedDiscreteTimeSystemMixin & operator=(const HypRedDiscreteTimeSystemMixin &) = delete;
  HypRedDiscreteTimeSystemMixin(HypRedDiscreteTimeSystemMixin &&) = default;
  HypRedDiscreteTimeSystemMixin & operator=(HypRedDiscreteTimeSystemMixin &&) = delete;
  ~HypRedDiscreteTimeSystemMixin() = default;

  template<class T1, class T2, class T3, class T4, class T5>
  HypRedDiscreteTimeSystemMixin(const T1 & romStateIn,
				const T2 & fomObj,
				const T3 & decoder,
				const T4 & fomNominalState,
				const T5 & projector)
    : T(romStateIn, fomObj, decoder, fomNominalState, projector),
      romSys_(romStateIn, T::projector_, T::fomCRef(), decoder, T::fomStatesMngr_)
  {}

  const rom_sys_t & romCRef() const{ return romSys_; }
};

template <class T, class rom_sys_t>
struct MaskedDiscreteTimeSystemMixin : T
{
  const rom_sys_t romSys_;

  MaskedDiscreteTimeSystemMixin() = delete;
  MaskedDiscreteTimeSystemMixin(const MaskedDiscreteTimeSystemMixin &) = default;
  MaskedDiscreteTimeSystemMixin & operator=(const MaskedDiscreteTimeSystemMixin &) = delete;
  MaskedDiscreteTimeSystemMixin(MaskedDiscreteTimeSystemMixin &&) = default;
  MaskedDiscreteTimeSystemMixin & operator=(MaskedDiscreteTimeSystemMixin &&) = delete;
  ~MaskedDiscreteTimeSystemMixin() = default;

  template<class T1, class T2, class T3, class T4, class T5, class T6>
  MaskedDiscreteTimeSystemMixin(const T1 & romStateIn,
				const T2 & fomObj,
				const T3 & decoder,
				const T4 & fomNominalState,
				const T5 & projector,
				const T6 & masker)
    : T(romStateIn, fomObj, decoder, fomNominalState, projector, masker),
      romSys_(romStateIn, T::projector_, T::masker_, T::fomCRef(), decoder, T::fomStatesMngr_)
  {}

  const rom_sys_t & romCRef() const{ return romSys_; }
};

template <class T, class rom_sys_t>
struct DefaultExplicitSystemMixin : T
{
  const rom_sys_t romSys_;

  DefaultExplicitSystemMixin() = delete;
  DefaultExplicitSystemMixin(const DefaultExplicitSystemMixin &) = default;
  DefaultExplicitSystemMixin & operator=(const DefaultExplicitSystemMixin &) = delete;
  DefaultExplicitSystemMixin(DefaultExplicitSystemMixin &&) = default;
  DefaultExplicitSystemMixin & operator=(DefaultExplicitSystemMixin &&) = delete;
  ~DefaultExplicitSystemMixin() = default;

  template<class T1, typename ...Args>
  DefaultExplicitSystemMixin(const T1 & romStateIn, Args && ...args)
    : T(romStateIn, std::forward<Args>(args)...),
      romSys_(romStateIn, T::projector_, T::fomCRef(), T::fomStatesMngr_)
  {}

  const rom_sys_t & romCRef() const{ return romSys_; }
};

template <class T, class rom_sys_t>
struct MaskedVeloExplicitSystemMixin : T
{
  const rom_sys_t romSys_;

  MaskedVeloExplicitSystemMixin() = delete;
  MaskedVeloExplicitSystemMixin(const MaskedVeloExplicitSystemMixin &) = default;
  MaskedVeloExplicitSystemMixin & operator=(const MaskedVeloExplicitSystemMixin &) = delete;
  MaskedVeloExplicitSystemMixin(MaskedVeloExplicitSystemMixin &&) = default;
  MaskedVeloExplicitSystemMixin & operator=(MaskedVeloExplicitSystemMixin &&) = delete;
  ~MaskedVeloExplicitSystemMixin() = default;

  template<class T1, typename ...Args>
  MaskedVeloExplicitSystemMixin(const T1 & romStateIn, Args && ...args)
    : T(romStateIn, std::forward<Args>(args)...),
      romSys_(romStateIn, T::projector_, T::masker_, T::fomCRef(), T::fomStatesMngr_)
  {}

  const rom_sys_t & romCRef() const{ return romSys_; }
};

template <class T, class rom_sys_t>
using HypRedVeloExplicitSystemMixin = DefaultExplicitSystemMixin<T, rom_sys_t>;

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
  ImplicitPoliciesMixin(const T1 & romStateIn,
			const T2 & fomObj,
			const T3 & decoder,
			const T4 & fomNominalState,
			Args && ... args)
    : T(romStateIn, fomObj, decoder, fomNominalState, std::forward<Args>(args)...),
      residualPolicy_(::pressio::ops::extent(romStateIn,0), T::projector_, T::fomCRef(), T::fomStatesMngr_),
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
  ImplicitPoliciesMixin(const T1 & romStateIn,
			const T2 & fomObj,
			const T3 & decoder,
			const T4 & fomNominalState,
			Args && ... args)
    : T(romStateIn, fomObj, decoder, fomNominalState, std::forward<Args>(args)...),
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
  ImplicitPoliciesMixin(const T1 & romStateIn,
			const T2 & fomObj,
			const T3 & decoder,
			const T4 & fomNominalState,
			Args && ... args)
    : T(romStateIn, fomObj, decoder, fomNominalState, std::forward<Args>(args)...),
      residualPolicy_(::pressio::ops::extent(romStateIn,0),
		      T::projector_, T::masker_, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(::pressio::ops::extent(romStateIn,0),
		      ::pressio::ops::extent(romStateIn,0),
		      T::projector_, T::masker_, T::fomCRef(), T::fomStatesMngr_, decoder)
  {}
};

// aliases to make things easier
template <class T, typename ...Args>
using DefaultImplicitPoliciesMixin = ImplicitPoliciesMixin<T, true, false, false, Args...>;

template <class T, typename ...Args>
using HypRedVeloImplicitPoliciesMixin = ImplicitPoliciesMixin<T, false, true, false, Args...>;

template <class T, typename ...Args>
using MaskedVeloImplicitPoliciesMixin = ImplicitPoliciesMixin<T, false, false, true, Args...>;

}}}}
#endif  // ROM_GALERKIN_IMPL_ROM_PROBLEM_MEMBERS_HPP_
