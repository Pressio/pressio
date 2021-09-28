/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_default_problem_explicit_stepping.hpp
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

#ifndef ROM_LSPG_IMPL_CONTINUOUS_TIME_API_ROM_LSPG_DEFAULT_PROBLEM_EXPLICIT_STEPPING_HPP_
#define ROM_LSPG_IMPL_CONTINUOUS_TIME_API_ROM_LSPG_DEFAULT_PROBLEM_EXPLICIT_STEPPING_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

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

template <class traits> struct SteadyMembersCommon
{
  using At = ::pressio::rom::impl::FomObjHolder<
    typename traits::fom_system_type>;

  using Bt = ::pressio::rom::lspg::impl::AddFomStatesManagerSteady<
    At,
    typename traits::fom_state_type,
    typename traits::fom_state_reconstr_type,
    typename traits::fom_states_manager_type>;
};

template <int flag, class traits>
struct SteadyMembers{ using type = void; };

// default
template <class traits> struct SteadyMembers<0, traits>
  : SteadyMembersCommon<traits>
{
  using base_t = SteadyMembersCommon<traits>;
  using typename base_t::Bt;

  using Ct = AddDefaultPolicies<
    Bt, typename traits::residual_policy_type, typename traits::jacobian_policy_type>;
  using type = SystemMixin<
    Ct, typename traits::steady_system_type>;
};

// masked
template <class traits> struct SteadyMembers<1, traits>
  : SteadyMembersCommon<traits>
{
  using base_t = SteadyMembersCommon<traits>;
  using typename base_t::Bt;

  using Ct = AddSinglyDecoratedPolicies<
    Bt,
    typename traits::masker_type,
    typename traits::residual_policy_type,
    typename traits::jacobian_policy_type>;
  using type = SystemMixin<
    Ct, typename traits::steady_system_type>;
};

// preconditioned default
template <class traits> struct SteadyMembers<2, traits>
  : SteadyMembersCommon<traits>
{
  using base_t = SteadyMembersCommon<traits>;
  using typename base_t::Bt;

  using Ct = AddSinglyDecoratedPolicies<
    Bt,
    typename traits::preconditioner_type,
    typename traits::residual_policy_type,
    typename traits::jacobian_policy_type>;
  using type = SystemMixin<
    Ct, typename traits::steady_system_type>;
};

// preconditioned masked
template <class traits> struct SteadyMembers<3, traits>
  : SteadyMembersCommon<traits>
{
  using base_t = SteadyMembersCommon<traits>;
  using typename base_t::Bt;

  using Ct = AddDoublyDecoratedPolicies<
    Bt,
    typename traits::preconditioner_type,
    typename traits::masker_type,
    typename traits::residual_policy_type,
    typename traits::jacobian_policy_type>;
  using type = SystemMixin<
    Ct, typename traits::steady_system_type>;
};

///////////////////
// problem class //
///////////////////
template <int flag, typename ...Args>
class SteadyProblem
{
public:
  using traits = ::pressio::Traits<SteadyProblem<flag, Args...>>;

private:
  using fom_system_type	 = typename traits::fom_system_type;
  using decoder_type	 = typename traits::decoder_type;
  using lspg_state_type = typename traits::lspg_state_type;
  using steady_system_type = typename traits::steady_system_type;
  using fom_state_type	 = typename traits::fom_state_type;
  using fom_state_reconstr_type = typename traits::fom_state_reconstr_type;
  typename SteadyMembers<flag, traits>::type members_;

public:
  steady_system_type & system(){ return members_.systemObj_; }

  const fom_state_type & currentFomState() const{
    return members_.fomStatesMngr_(::pressio::ode::n());
  }

  const fom_state_reconstr_type & fomStateReconstructor() const{
    return members_.fomStateReconstructor_;
  }

public:
  SteadyProblem() = delete;
  SteadyProblem(const SteadyProblem &) = default;
  SteadyProblem & operator=(const SteadyProblem &) = delete;
  SteadyProblem(SteadyProblem &&) = default;
  SteadyProblem & operator=(SteadyProblem &&) = delete;
  ~SteadyProblem() = default;

  template<
    int _flag = flag,
    mpl::enable_if_t<_flag==0, int> = 0
    >
  SteadyProblem(const fom_system_type & fomObj,
		decoder_type & decoder,
		const lspg_state_type & romState,
		const fom_state_type & fomNominalState)
    : members_(romState, fomObj, decoder, fomNominalState){}

  template<
    int _flag = flag, class ...Args2,
    mpl::enable_if_t<_flag==1 or _flag==2 or _flag==3, int> = 0
    >
  SteadyProblem(const fom_system_type & fomObj,
		decoder_type & decoder,
		const lspg_state_type & romState,
		const fom_state_type & fomNominalState,
		Args2 && ...args)
    : members_(romState, fomObj, decoder,
	       fomNominalState, std::forward<Args2>(args) ...){}
};

}}}}//end namespace pressio::rom::lspg::impl
#endif  // ROM_LSPG_IMPL_CONTINUOUS_TIME_API_ROM_LSPG_DEFAULT_PROBLEM_EXPLICIT_STEPPING_HPP_
