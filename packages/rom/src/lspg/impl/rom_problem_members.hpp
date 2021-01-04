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

#ifndef ROM_LSPG_IMPL_ROM_PROBLEM_MEMBERS_HPP_
#define ROM_LSPG_IMPL_ROM_PROBLEM_MEMBERS_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template <
  typename T,
  bool def, bool masked, bool prec, bool hypred,
  typename ... Args
  >
struct PoliciesMixin;

// specialize for default
template <class T, class ops_t, class r_pol_t, class j_pol_t>
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
    class T1, class T2, class T3, class T4, class _ops_t = ops_t,
    mpl::enable_if_t<std::is_void<_ops_t>::value, int > = 0
    >
  PoliciesMixin(const T1 & romStateIn,
		const T2 & fomObj,
		T3 & decoder,
		const T4 & fomNominalStateNative)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative),
      residualPolicy_(T::fomStatesMngr_),
      jacobianPolicy_(T::fomStatesMngr_, decoder)
  {}

  template<
    class T1, class T2, class T3, class T4, class _ops_t = ops_t,
    mpl::enable_if_t<!std::is_void<_ops_t>::value, int > = 0
    >
  PoliciesMixin(const T1 & romStateIn,
		const T2 & fomObj,
		T3 & decoder,
		const T4 & fomNominalStateNative,
		const _ops_t & ops)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative, ops),
      residualPolicy_(T::fomStatesMngr_, ops),
      jacobianPolicy_(T::fomStatesMngr_, decoder, ops)
  {}
};

// specialize for masked
template <class T, class masker_t, class ops_t, class r_pol_t, class j_pol_t>
struct PoliciesMixin<
  T, false, true, false, false, masker_t, ops_t, r_pol_t, j_pol_t
  > : T
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

  PoliciesMixin() = delete;
  PoliciesMixin(const PoliciesMixin &) = default;
  PoliciesMixin & operator=(const PoliciesMixin &) = delete;
  PoliciesMixin(PoliciesMixin &&) = default;
  PoliciesMixin & operator=(PoliciesMixin &&) = delete;
  ~PoliciesMixin() = default;

  template<
    class T1, class T2, class T3, class T4, class T5,
    class _ops_t = ops_t,
    mpl::enable_if_t<std::is_void<_ops_t>::value, int > = 0
    >
  PoliciesMixin(const T1 & romStateIn,
		const T2 & fomObj,
		T3 & decoder,
		const T4 & fomNominalStateNative,
		const T5 & masker)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative),
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
      masker_(masker),
      residualPolicy_(masker_, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(masker_, T::fomCRef(), T::fomStatesMngr_, decoder)
#else
      residualPolicy_(masker, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(masker, T::fomCRef(), T::fomStatesMngr_, decoder)
#endif
  {
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
    static_assert
      (std::is_same<T5, pybind11::object>::value,
       "Maked policies mixin: masker object must be a pybind11::object");
#endif
  }
};

// specialize for precond default
template <class T, class ops_t, class r_pol_t, class j_pol_t>
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
    class T1, class T2, class T3, class T4, class T5,
    class _ops_t = ops_t,
    mpl::enable_if_t<std::is_void<_ops_t>::value, int > = 0
    >
  PoliciesMixin(const T1 & romStateIn,
		const T2 & fomObj,
		T3 & decoder,
		const T4 & fomNominalStateNative,
		const T5 & preconditioner)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative),
      residualPolicy_(preconditioner, T::fomStatesMngr_),
      jacobianPolicy_(preconditioner, T::fomStatesMngr_, decoder)
  {}
};

// specialize for hyp-red with void stencil-to-sample mapping
template <class T, class ops_t, class r_pol_t, class j_pol_t>
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
    class T1, class T2, class T3, class T4,
    class _ops_t = ops_t,
    mpl::enable_if_t<std::is_void<_ops_t>::value, int > = 0
    >
  PoliciesMixin(const T1 & romStateIn,
		const T2 & fomObj,
		T3 & decoder,
		const T4 & fomNominalStateNative)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative),
      residualPolicy_(T::fomStatesMngr_),
      jacobianPolicy_(T::fomStatesMngr_, decoder)
  {}
};

// specialize for hyp-red with nonvoid stencil-to-sample mapping
template <class T, class ops_t, class r_pol_t, class j_pol_t, class sTos_t>
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
    class T1, class T2, class T3, class T4, class T5,
    class _ops_t = ops_t,
    mpl::enable_if_t<std::is_void<_ops_t>::value, int > = 0
    >
  PoliciesMixin(const T1 & romStateIn,
		const T2 & fomObj,
		T3 & decoder,
		const T4 & fomNominalStateNative,
		const T5 & meshToStencilMapper)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative),
      meshToStencilMapper_(meshToStencilMapper),
      residualPolicy_(T::fomStatesMngr_, meshToStencilMapper_),
      jacobianPolicy_(T::fomStatesMngr_, decoder, meshToStencilMapper_)
  {}
};

// specialize for preconditioned hyp-red with void stencil-to-sample mapping
template <class T, class ops_t, class r_pol_t, class j_pol_t>
struct PoliciesMixin<
  T, false, false, true, true, ops_t, r_pol_t, j_pol_t, void
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
    class T1, class T2, class T3, class T4, class T5,
    class _ops_t = ops_t,
    mpl::enable_if_t<std::is_void<_ops_t>::value, int > = 0
    >
  PoliciesMixin(const T1 & romStateIn,
		const T2 & fomObj,
		T3 & decoder,
		const T4 & fomNominalStateNative,
		const T5 & preconditioner)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative),
      residualPolicy_(preconditioner, T::fomStatesMngr_),
      jacobianPolicy_(preconditioner, T::fomStatesMngr_, decoder)
  {}
};

// specialize for preconditioned hyp-red with nonvoid stencil-to-sample mapping
template <class T, class ops_t, class r_pol_t, class j_pol_t, class sTos_t>
struct PoliciesMixin<
  T, false, false, true, true, ops_t, r_pol_t, j_pol_t, sTos_t
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
    class T1, class T2, class T3, class T4, class T5, class T6,
    class _ops_t = ops_t,
    mpl::enable_if_t<std::is_void<_ops_t>::value, int > = 0
    >
  PoliciesMixin(const T1 & romStateIn,
		const T2 & fomObj,
		T3 & decoder,
		const T4 & fomNominalStateNative,
		const T5 & preconditioner,
		const T6 & meshToStencilMapper)
    : T(fomObj, decoder, romStateIn, fomNominalStateNative),
      meshToStencilMapper_(meshToStencilMapper),
      residualPolicy_(preconditioner, T::fomStatesMngr_, meshToStencilMapper_),
      jacobianPolicy_(preconditioner, T::fomStatesMngr_, decoder, meshToStencilMapper_)
  {}
};

// aliases to make things easier
template <class T, typename ...Args>
using DefaultPoliciesMixin = PoliciesMixin<T, true, false, false, false, Args...>;

template <class T, typename ...Args>
using MaskedPoliciesMixin = PoliciesMixin<T, false, true, false, false, Args...>;

template <class T, typename ...Args>
using PrecondPoliciesMixin = PoliciesMixin<T, false, false, true, false, Args...>;

template <class T, typename ...Args>
using HypRedPoliciesMixin = PoliciesMixin<T, false, false, false, true, Args...>;

template <class T, typename ...Args>
using PrecHypRedPoliciesMixin = PoliciesMixin<T, false, false, true, true, Args...>;

//---------------------------------------------------
// system mixin (used for steady LSPG)
//---------------------------------------------------
template <class T, class system_t>
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
#endif  // ROM_LSPG_IMPL_ROM_PROBLEM_MEMBERS_HPP_
