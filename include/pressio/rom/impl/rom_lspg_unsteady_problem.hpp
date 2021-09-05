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

#ifndef ROM_LSPG_IMPL_UNSTEADY_PROBLEM_HPP_
#define ROM_LSPG_IMPL_UNSTEADY_PROBLEM_HPP_

namespace pressio{ namespace rom{ namespace lspg{ namespace impl{

template <class traits> struct UnsteadyMembersCommon
{
  static constexpr auto binding_sentinel = traits::binding_sentinel;
  using At = ::pressio::rom::impl::FomObjMixin<
    typename traits::fom_system_type, binding_sentinel>;

  using Bt = ::pressio::rom::impl::FomStatesMngrMixin<
    At,
    typename traits::fom_state_type,
    typename traits::fom_state_reconstr_type,
    typename traits::fom_states_manager_type>;
};

template <int flag, class traits>
struct UnsteadyMembers{ using type = void; };

// default
template <class traits> struct UnsteadyMembers<0, traits> : UnsteadyMembersCommon<traits>
{
  using base_t = UnsteadyMembersCommon<traits>;
  using typename base_t::At;
  using typename base_t::Bt;
  using base_t::binding_sentinel;

  using Ct  = DefaultPoliciesMixin<
    Bt,
    typename traits::residual_policy_type,
    typename traits::jacobian_policy_type>;
  using type = ::pressio::rom::impl::ImplicitStepperMixin<
    Ct, typename traits::stepper_type>;
};

///////////////////
// problem class //
///////////////////
template <int flag, typename ...Args>
class UnsteadyProblem
{
public:
  using traits = ::pressio::Traits<UnsteadyProblem<flag, Args...>>;

private:
  using fom_system_type	 = typename traits::fom_system_type;
  using decoder_type	 = typename traits::decoder_type;
  using lspg_state_type = typename traits::lspg_state_type;
  using stepper_type = typename traits::stepper_type;
  using fom_state_type	 = typename traits::fom_state_type;
  using fom_state_reconstr_type = typename traits::fom_state_reconstr_type;
  typename UnsteadyMembers<flag, traits>::type members_;

public:
  stepper_type & stepper(){ return members_.stepperObj_; }

  const fom_state_type & currentFomState() const{
    return members_.fomStatesMngr_(::pressio::ode::n());
  }

  const fom_state_reconstr_type & fomStateReconstructor() const{
    return members_.fomStateReconstructor_;
  }

public:
  UnsteadyProblem() = delete;
  UnsteadyProblem(const UnsteadyProblem &) = default;
  UnsteadyProblem & operator=(const UnsteadyProblem &) = delete;
  UnsteadyProblem(UnsteadyProblem &&) = default;
  UnsteadyProblem & operator=(UnsteadyProblem &&) = delete;
  ~UnsteadyProblem() = default;

  template<
    int _flag = flag,
    mpl::enable_if_t<_flag==0, int> = 0
    >
  UnsteadyProblem(const fom_system_type & fomObj,
		decoder_type & decoder,
		const lspg_state_type & romState,
		const fom_state_type & fomNominalState)
    : members_(romState, fomObj, decoder, fomNominalState){}
};

}}}}//end namespace pressio::rom::lspg::impl
#endif  // ROM_LSPG_IMPL_CONTINUOUS_TIME_API_ROM_LSPG_DEFAULT_PROBLEM_EXPLICIT_STEPPING_HPP_
