/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_problem.hpp
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

#ifndef ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_DEFAULT_PROBLEM_EXPLICIT_STEPPING_HPP_
#define ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_DEFAULT_PROBLEM_EXPLICIT_STEPPING_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <class traits> struct MembersCommon
{
  using At = ::pressio::rom::impl::FomObjHolder<
    typename traits::fom_system_type, traits::binding_sentinel>;

  using Bt = ::pressio::rom::galerkin::impl::AddFomStatesManager<
    At,
    traits::is_cont_time,
    typename traits::fom_state_type,
    typename traits::fom_state_reconstr_type,
    typename traits::fom_states_manager_type>;
};

template <int flag, class traits>
struct Members{ using type = void; };

// --------
// default
// --------
// cont-time explicit
template <class traits> struct Members<0, traits> : MembersCommon<traits>
{
  using base_t = MembersCommon<traits>;
  using typename base_t::Bt;

  using Ct = AddProjector<Bt, typename traits::projector_type>;
  using Dt = AddDefaultExplicitSystem<Ct, typename traits::rom_system_type>;
  using type = ::pressio::rom::impl::AddExplicitStepper<
    Dt, typename traits::stepper_type>;
};

// cont-time implicit
template <class traits> struct Members<1, traits> : MembersCommon<traits>
{
  using base_t = MembersCommon<traits>;
  using typename base_t::Bt;

  using Ct  = AddProjector<Bt, typename traits::projector_type>;
  using Dt  = AddDefaultImplicitPolicies<
    Ct,
    typename traits::residual_policy_type,
    typename traits::jacobian_policy_type>;
  using type = ::pressio::rom::impl::AddImplicitStepper<
    Dt, typename traits::stepper_type>;
};

// disc-time implicit
template <class traits> struct Members<2, traits> : MembersCommon<traits>
{
  using base_t = MembersCommon<traits>;
  using typename base_t::Bt;

  using Ct  = AddProjector<Bt, typename traits::projector_type>;
  using Dt  = AddDefaultDiscreteTimeSystem<Ct, typename traits::rom_system_type>;
  using type = ::pressio::rom::impl::AddImplicitArbStepper<
    Dt, typename traits::stepper_type>;
};

// --------
// masked
// --------
// cont-time explicit
template <class traits> struct Members<3, traits> : MembersCommon<traits>
{
  using base_t = MembersCommon<traits>;
  using typename base_t::Bt;

  using Ct   = AddMasker<Bt, typename traits::masker_type>;
  using Dt   = AddProjector<Ct, typename traits::projector_type>;
  using Et   = AddMaskedVeloExplicitSystem<Dt, typename traits::rom_system_type>;
  using type = ::pressio::rom::impl::AddExplicitStepper<
    Et, typename traits::stepper_type>;
};

// cont-time implicit
template <class traits> struct Members<4, traits> : MembersCommon<traits>
{
  using base_t = MembersCommon<traits>;
  using typename base_t::Bt;

  using Ct  = AddMasker<Bt, typename traits::masker_type>;
  using Dt  = AddProjector<Ct, typename traits::projector_type>;
  using Et  = AddMaskedVeloImplicitPolicies<
    Dt, typename traits::residual_policy_type, typename traits::jacobian_policy_type>;
  using type= ::pressio::rom::impl::AddImplicitStepper<
    Et, typename traits::stepper_type>;
};

// disc-time masked
template <class traits> struct Members<5, traits> : MembersCommon<traits>
{
  using base_t = MembersCommon<traits>;
  using typename base_t::Bt;

  using Ct  = AddMasker<Bt, typename traits::masker_type>;
  using Dt  = AddProjector<Ct, typename traits::projector_type>;
  using Et  = AddMaskedDiscreteTimeSystem<Dt, typename traits::rom_system_type>;
  using type = ::pressio::rom::impl::AddImplicitArbStepper<
    Et, typename traits::stepper_type>;
};

// --------
// hyp-red
// --------
// cont-time explicit
template <class traits> struct Members<6, traits> : MembersCommon<traits>
{
  using base_t = MembersCommon<traits>;
  using typename base_t::Bt;

  using Ct = AddProjector<Bt, typename traits::projector_type>;
  using Dt = AddHypRedVeloExplicitSystem<Ct, typename traits::rom_system_type>;
  using type = ::pressio::rom::impl::AddExplicitStepper<
    Dt, typename traits::stepper_type>;
};

// cont-time implicit
template <class traits> struct Members<7, traits> : MembersCommon<traits>
{
  using base_t = MembersCommon<traits>;
  using typename base_t::Bt;

  using Ct  = AddProjector<Bt, typename traits::projector_type>;
  using Dt  = AddHypRedVeloImplicitPolicies<
    Ct, typename traits::residual_policy_type, typename traits::jacobian_policy_type>;
  using type= ::pressio::rom::impl::AddImplicitStepper<
    Dt, typename traits::stepper_type>;
};

// disc-time implicit
template <class traits> struct Members<8, traits> : MembersCommon<traits>
{
  using base_t = MembersCommon<traits>;
  using typename base_t::Bt;

  using Ct = AddProjector<Bt, typename traits::projector_type>;
  using Dt = AddHypRedDiscreteTimeSystem<Ct, typename traits::rom_system_type>;
  using type = ::pressio::rom::impl::AddImplicitArbStepper<
    Dt, typename traits::stepper_type>;
};

///////////////////
// problem class //
///////////////////
template <int flag, typename ...Args>
class Problem
{
public:
  using traits = ::pressio::Traits<Problem<flag, Args...>>;

private:
  using fom_system_type	 = typename traits::fom_system_type;
  using decoder_type	 = typename traits::decoder_type;
  using galerkin_state_type = typename traits::galerkin_state_type;
  using stepper_type     = typename traits::stepper_type;
  using fom_state_type	 = typename traits::fom_state_type;
  using fom_state_reconstr_type = typename traits::fom_state_reconstr_type;
  typename Members<flag, traits>::type members_;

public:
  stepper_type & stepper(){ return members_.stepperObj_; }

  const fom_state_type & currentFomState() const{
    return members_.fomStatesMngr_(::pressio::ode::n());
  }

  const fom_state_reconstr_type & fomStateReconstructor() const{
    return members_.fomStateReconstructor_;
  }

public:
  Problem() = delete;
  Problem(const Problem &) = default;
  Problem & operator=(const Problem &) = delete;
  Problem(Problem &&) = default;
  Problem & operator=(Problem &&) = delete;
  ~Problem() = default;

  template<
    int _flag = flag,
    mpl::enable_if_t<_flag<=2, int> = 0
    >
  Problem(::pressio::ode::StepScheme name,
	  const fom_system_type & fomObj,
	  decoder_type & decoder,
	  const galerkin_state_type & romState,
	  const fom_state_type & fomNominalState)
    : members_(name, romState, fomObj, decoder, fomNominalState){}

  template<
    int _flag = flag, class ...Args2,
    mpl::enable_if_t<_flag>=3 and _flag<=5, int> = 0
    >
  Problem(::pressio::ode::StepScheme name,
	  const fom_system_type & fomObj,
	  decoder_type & decoder,
	  const galerkin_state_type & romState,
	  const fom_state_type & fomNominalState,
	  const typename traits::projector_type & projector,
	  Args2 && ...args)
    : members_(name, romState, fomObj, decoder, fomNominalState,
	       projector, std::forward<Args2>(args) ...){}

  template<
    int _flag = flag, class ...Args2,
    mpl::enable_if_t<_flag>=6 and _flag<=8, int> = 0
    >
  Problem(::pressio::ode::StepScheme name,
	  const fom_system_type & fomObj,
	  decoder_type & decoder,
	  const galerkin_state_type & romState,
	  const fom_state_type & fomNominalState,
	  const typename traits::projector_type & projector)
    : members_(name, romState, fomObj, decoder, fomNominalState, projector){}
};

}}}}//end namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_DEFAULT_PROBLEM_EXPLICIT_STEPPING_HPP_
