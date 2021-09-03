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

// template <
//   class T, bool def, bool masked, bool prec, bool hypred,
//   class ... Args
//   >
// struct PoliciesMixin;

template <class T, class r_pol_t, class j_pol_t>
struct DefaultPoliciesMixin : T
{
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  DefaultPoliciesMixin() = delete;
  DefaultPoliciesMixin(const DefaultPoliciesMixin &) = default;
  DefaultPoliciesMixin & operator=(const DefaultPoliciesMixin &) = delete;
  DefaultPoliciesMixin(DefaultPoliciesMixin &&) = default;
  DefaultPoliciesMixin & operator=(DefaultPoliciesMixin &&) = delete;
  ~DefaultPoliciesMixin() = default;

  template<class T1, class T2, class T3, class T4>
  DefaultPoliciesMixin(const T1 & romStateIn,
		       const T2 & fomObj,
		       T3 & decoder,
		       const T4 & fomNominalState)
    : T(fomObj, decoder, romStateIn, fomNominalState),
      residualPolicy_(T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(T::fomCRef(), T::fomStatesMngr_, decoder)
  {}
};

template <class T, class UserProvidedFunctor_t, class r_pol_t, class j_pol_t>
struct SinglyDecoratedPoliciesMixin : T
{
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  SinglyDecoratedPoliciesMixin() = delete;
  SinglyDecoratedPoliciesMixin(const SinglyDecoratedPoliciesMixin &) = default;
  SinglyDecoratedPoliciesMixin & operator=(const SinglyDecoratedPoliciesMixin &) = delete;
  SinglyDecoratedPoliciesMixin(SinglyDecoratedPoliciesMixin &&) = default;
  SinglyDecoratedPoliciesMixin & operator=(SinglyDecoratedPoliciesMixin &&) = delete;
  ~SinglyDecoratedPoliciesMixin() = default;

  template<class T1, class T2, class T3, class T4>
  SinglyDecoratedPoliciesMixin(const T1 & romStateIn,
			       const T2 & fomObj,
			       T3 & decoder,
			       const T4 & fomNominalState,
			       const UserProvidedFunctor_t & userProvidedFunctor)
    : T(fomObj, decoder, romStateIn, fomNominalState),
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
struct DoublyDecoratedPoliciesMixin : T
{
  r_pol_t residualPolicy_;
  j_pol_t jacobianPolicy_;

  DoublyDecoratedPoliciesMixin() = delete;
  DoublyDecoratedPoliciesMixin(const DoublyDecoratedPoliciesMixin &) = default;
  DoublyDecoratedPoliciesMixin & operator=(const DoublyDecoratedPoliciesMixin &) = delete;
  DoublyDecoratedPoliciesMixin(DoublyDecoratedPoliciesMixin &&) = default;
  DoublyDecoratedPoliciesMixin & operator=(DoublyDecoratedPoliciesMixin &&) = delete;
  ~DoublyDecoratedPoliciesMixin() = default;

  template<class T1, class T2, class T3, class T4>
  DoublyDecoratedPoliciesMixin(const T1 & romStateIn,
			       const T2 & fomObj,
			       T3 & decoder,
			       const T4 & fomNominalState,
			       const UserProvidedFunctor1_t & functor1,
			       const UserProvidedFunctor2_t & functor2)
    : T(fomObj, decoder, romStateIn, fomNominalState),
      residualPolicy_(functor1, functor2, T::fomCRef(), T::fomStatesMngr_),
      jacobianPolicy_(functor1, functor2, T::fomCRef(), T::fomStatesMngr_, decoder)
  {}
};

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

  template<class...Args>
  SystemMixin(Args && ...args)
    : T(std::forward<Args>(args)...),
      systemObj_(T::residualPolicy_, T::jacobianPolicy_)
  {}
};

}}}}
#endif  // ROM_LSPG_IMPL_ROM_PROBLEM_MEMBERS_HPP_
