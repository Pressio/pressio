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

template <class T, class ops_t, class projector_t>
struct ProjectorMixin : T
{
  projector_t projector_;

  ProjectorMixin() = delete;
  ProjectorMixin(const ProjectorMixin &) = default;
  ProjectorMixin & operator=(const ProjectorMixin &) = delete;
  ProjectorMixin(ProjectorMixin &&) = default;
  ProjectorMixin & operator=(ProjectorMixin &&) = delete;
  ~ProjectorMixin() = default;

  template<
    class T1, class T2, class T3, class T4, class _ops_t = ops_t,
    mpl::enable_if_t<std::is_void<_ops_t>::value, int > = 0
    >
  ProjectorMixin(const T1 & romStateIn,
		 const T2 & fomObj,
		 const T3 & decoder,
		 const T4 & fomNominalStateNative)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative),
      projector_(decoder)
  {}

  template<
    class T1, class T2, class T3, class T4, class _ops_t = ops_t,
    mpl::enable_if_t<mpl::not_void<_ops_t>::value, int > = 0
    >
  ProjectorMixin(const T1 & romStateIn,
		 const T2 & fomObj,
		 const T3 & decoder,
		 const T4 & fomNominalStateNative,
		 const _ops_t & udOps)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative, udOps),
      projector_(decoder, udOps)
  {}
};

//---------------------------------------------------
// implicit policies
//---------------------------------------------------
template <class T, bool def, bool hypredvelo, bool, typename ... Args>
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
			const T4 & fomNominalStateNative,
			Args && ... args)
    : T(romStateIn, fomObj, decoder, fomNominalStateNative, std::forward<Args>(args)...),
      residualPolicy_(romStateIn.extent(0), T::projector_, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(romStateIn.extent(0), T::projector_, T::fomCRef(), T::fomStatesMngr_, decoder)
  {}
};

// specialize for hyp-red velo
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

  template<class T1, class T2, class T3, class T4, class T5, typename ...Args>
  ImplicitPoliciesMixin(const T1 & romStateIn,
			const T2 & fomObj,
			const T3 & decoder,
			const T4 & fomNominalStateNative,
			const T5 & projector,
			Args && ... args)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative, std::forward<Args>(args)...),
      residualPolicy_(romStateIn.extent(0), projector, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(romStateIn.extent(0), projector, T::fomCRef(), T::fomStatesMngr_, decoder)
  {}
};

// specialize for masked velo
template <class T, class masker_t, class r_pol_t, class j_pol_t>
struct ImplicitPoliciesMixin<T, false, false, true, masker_t, r_pol_t, j_pol_t> : T
{
/* here we need to consider also the case where the masker is a pybind11:object
   that is passed in directly from python: in that scenario, masker_t is
   a C++ wrapper class wrapping the actualy pure python class,
   so we need to create an object of this masker_t and pass that to the policies
   because the policies do NOT accept pybind11::objects
*/
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  const masker_t masker_;
#endif
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  ImplicitPoliciesMixin() = delete;
  ImplicitPoliciesMixin(const ImplicitPoliciesMixin &) = default;
  ImplicitPoliciesMixin & operator=(const ImplicitPoliciesMixin &) = delete;
  ImplicitPoliciesMixin(ImplicitPoliciesMixin &&) = default;
  ImplicitPoliciesMixin & operator=(ImplicitPoliciesMixin &&) = delete;
  ~ImplicitPoliciesMixin() = default;

  template<class T1, class T2, class T3, class T4, class T5, class T6, class ...Args>
  ImplicitPoliciesMixin(const T1 & romStateIn,
			const T2 & fomObj,
			const T3 & decoder,
			const T4 & fomNominalStateNative,
			const T5 & masker,
			const T6 & projector,
			Args && ... args)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative, std::forward<Args>(args)...),
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
      masker_(masker),
      residualPolicy_(romStateIn.extent(0), projector, masker_, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(romStateIn.extent(0), projector, masker_, T::fomCRef(), T::fomStatesMngr_, decoder)
#else
      residualPolicy_(romStateIn.extent(0), projector, masker, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(romStateIn.extent(0), projector, masker, T::fomCRef(), T::fomStatesMngr_, decoder)
#endif
  {}
};

// aliases to make things easier
template <class T, typename ...Args>
using DefaultImplicitPoliciesMixin = ImplicitPoliciesMixin<T, true, false, false, Args...>;

template <class T, typename ...Args>
using HypRedVeloImplicitPoliciesMixin = ImplicitPoliciesMixin<T, false, true, false, Args...>;
template <class T, typename ...Args>
using HypRedResidualImplicitPoliciesMixin = ImplicitPoliciesMixin<T, false, true, false, Args...>;

template <class T, typename ...Args>
using MaskedVeloImplicitPoliciesMixin = ImplicitPoliciesMixin<T, false, false, true, Args...>;
template <class T, typename ...Args>
using MaskedResidualImplicitPoliciesMixin = ImplicitPoliciesMixin<T, false, false, true, Args...>;


//---------------------------------------------------
// explicit policies
//---------------------------------------------------
template <class T, bool, bool, bool, typename ... Args>
struct ExplicitPoliciesMixin;

// default
template <class T, class rhs_pol_t>
struct ExplicitPoliciesMixin<T, true, false, false, rhs_pol_t> : T
{
  rhs_pol_t rhsPolicy_;

  ExplicitPoliciesMixin() = delete;
  ExplicitPoliciesMixin(const ExplicitPoliciesMixin &) = default;
  ExplicitPoliciesMixin & operator=(const ExplicitPoliciesMixin &) = delete;
  ExplicitPoliciesMixin(ExplicitPoliciesMixin &&) = default;
  ExplicitPoliciesMixin & operator=(ExplicitPoliciesMixin &&) = delete;
  ~ExplicitPoliciesMixin() = default;

  template<
    class T1, typename ...Args,
    mpl::enable_if_t<::pressio::containers::details::traits<T1>::rank==1,int> = 0
    >
  ExplicitPoliciesMixin(const T1 & romStateIn, Args && ...args)
    : T(romStateIn, std::forward<Args>(args)...),
      rhsPolicy_(romStateIn.extent(0), T::projector_, T::fomCRef(), T::fomStatesMngr_)
  {}

  template<
    class T1, typename ...Args,
    mpl::enable_if_t<::pressio::containers::details::traits<T1>::rank==2,int> = 0
    >
  ExplicitPoliciesMixin(const T1 & romStateIn, Args && ...args)
    : T(romStateIn, std::forward<Args>(args)...),
      rhsPolicy_(romStateIn.extent(0), romStateIn.extent(1),
		 T::projector_, T::fomCRef(), T::fomStatesMngr_)
  {}
};

// hypere-reduced
template <class T, class rhs_pol_t>
struct ExplicitPoliciesMixin<T, false, true, false, rhs_pol_t> : T
{
  rhs_pol_t rhsPolicy_;

  ExplicitPoliciesMixin() = delete;
  ExplicitPoliciesMixin(const ExplicitPoliciesMixin &) = default;
  ExplicitPoliciesMixin & operator=(const ExplicitPoliciesMixin &) = delete;
  ExplicitPoliciesMixin(ExplicitPoliciesMixin &&) = default;
  ExplicitPoliciesMixin & operator=(ExplicitPoliciesMixin &&) = delete;
  ~ExplicitPoliciesMixin() = default;

  template<
    class T1, class T2, class T3, class T4, class T5, typename ...Args,
    mpl::enable_if_t<::pressio::containers::details::traits<T1>::rank==1,int> = 0
    >
  ExplicitPoliciesMixin(const T1 & romStateIn,
			const T2 & fomObj,
			const T3 & decoder,
			const T4 & fomNominalStateNative,
			const T5 & projector,
			Args && ...args)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative, std::forward<Args>(args)...),
      rhsPolicy_(romStateIn.extent(0), projector, T::fomCRef(), T::fomStatesMngr_)
  {}

  template<
    class T1, class T2, class T3, class T4, class T5, typename ...Args,
    mpl::enable_if_t<::pressio::containers::details::traits<T1>::rank==2,int> = 0
    >
  ExplicitPoliciesMixin(const T1 & romStateIn,
			const T2 & fomObj,
			const T3 & decoder,
			const T4 & fomNominalStateNative,
			const T5 & projector,
			Args && ...args)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative, std::forward<Args>(args)...),
      rhsPolicy_(romStateIn.extent(0), romStateIn.extent(1),
		 projector, T::fomCRef(), T::fomStatesMngr_)
  {}
};

// masked
template <class T, class masker_t, class rhs_pol_t>
struct ExplicitPoliciesMixin<T, false, false, true, masker_t, rhs_pol_t> : T
{
/* here we need to consider also the case where the masker is a pybind11:object
   that is passed in directly from python: in that scenario, masker_t is
   a C++ wrapper class wrapping the actualy pure python class,
   so we need to create an object of this masker_t and pass that to the policies
   because the policies do NOT accept pybind11::objects
*/
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  const masker_t masker_;
#endif
  rhs_pol_t rhsPolicy_;

  ExplicitPoliciesMixin() = delete;
  ExplicitPoliciesMixin(const ExplicitPoliciesMixin &) = default;
  ExplicitPoliciesMixin & operator=(const ExplicitPoliciesMixin &) = delete;
  ExplicitPoliciesMixin(ExplicitPoliciesMixin &&) = default;
  ExplicitPoliciesMixin & operator=(ExplicitPoliciesMixin &&) = delete;
  ~ExplicitPoliciesMixin() = default;

  template<
    class T1, class T2, class T3, class T4, class T5, class T6, class ...Args,
    mpl::enable_if_t<::pressio::containers::details::traits<T1>::rank==1,int> = 0
    >
  ExplicitPoliciesMixin(const T1 & romStateIn,
			const T2 & fomObj,
			const T3 & decoder,
			const T4 & fomNominalStateNative,
			const T5 & masker,
			const T6 & projector,
			Args && ...args)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative, std::forward<Args>(args)...),
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
      masker_(masker),
      rhsPolicy_(romStateIn.extent(0), projector, masker_, T::fomCRef(), T::fomStatesMngr_)
#else
      rhsPolicy_(romStateIn.extent(0), projector, masker, T::fomCRef(), T::fomStatesMngr_)
#endif
  {
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    static_assert
      (std::is_same<T5, pybind11::object>::value,
       "Maked policies mixin: masker object must be a pybind11::object");
#endif
  }

  template<
    class T1, class T2, class T3, class T4, class T5, class T6, class ...Args,
    mpl::enable_if_t<::pressio::containers::details::traits<T1>::rank==2,int> = 0
    >
  ExplicitPoliciesMixin(const T1 & romStateIn,
			const T2 & fomObj,
			const T3 & decoder,
			const T4 & fomNominalStateNative,
			const T5 & masker,
			const T6 & projector,
			Args && ...args)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative, std::forward<Args>(args)...),
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
      masker_(masker),
      rhsPolicy_(romStateIn.extent(0), romStateIn.extent(1),
		 projector, masker_, T::fomCRef(), T::fomStatesMngr_)
#else
      rhsPolicy_(romStateIn.extent(0), romStateIn.extent(1),
		 projector, masker, T::fomCRef(), T::fomStatesMngr_)
#endif
  {
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    static_assert
      (std::is_same<T5, pybind11::object>::value,
       "Maked policies mixin: masker object must be a pybind11::object");
#endif
  }
};

// aliases to make things easier
template <class T, typename ...Args>
using DefaultExplicitPoliciesMixin = ExplicitPoliciesMixin<T, true, false, false, Args...>;

template <class T, typename ...Args>
using HypRedVeloExplicitPoliciesMixin = ExplicitPoliciesMixin<T, false, true, false, Args...>;

template <class T, typename ...Args>
using MaskedVeloExplicitPoliciesMixin = ExplicitPoliciesMixin<T, false, false, true, Args...>;

}}}}
#endif  // ROM_GALERKIN_IMPL_ROM_PROBLEM_MEMBERS_HPP_
