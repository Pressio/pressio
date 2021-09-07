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
  class T,
  bool is_static,
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

  // template<
  //   class T1, class T2, class T3, class T4,
  //   bool _is_static = is_static,
  //   mpl::enable_if_t<!_is_static, int> = 0
  //   >
  // AddFomStatesManager(::pressio::ode::SteppersE name,
  // 		      const T1 & fomObj,
  // 		      const T2 & decoder,
  // 		      const T3 & romStateIn,
  // 		      const T4 & fomNominalStateNative)
  //   : T(fomObj),
  //     fomNominalState_(::pressio::ops::clone(fomNominalStateNative)),
  //     fomStateReconstructor_(fomNominalState_, decoder),
  //     fomStatesMngr_(create_manager_stencil_fom_states<
  // 		     ManagerStencilFomStatesType>(name, fomStateReconstructor_, fomNominalState_)
  // 		     )
  // {
  //   // reconstruct current fom state so that we have something
  //   // consisten with the current romState
  //   fomStatesMngr_.reconstructCurrentFomState(romStateIn);
  // }

  template<
    class T1, class T2, class T3, class T4,
    bool _is_static = is_static,
    mpl::enable_if_t<_is_static, int> = 0
    >
  AddFomStatesManager(const T1 & fomObj,
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
};

}}}}
#endif  // ROM_LSPG_IMPL_ROM_PROBLEM_MEMBERS_HPP_
