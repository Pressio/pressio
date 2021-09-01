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

#ifndef ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_DEFAULT_PROBLEM_EXPLICIT_STEPPING_HPP_
#define ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_DEFAULT_PROBLEM_EXPLICIT_STEPPING_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <int flag, class traits>
struct Members{ using type = void; };

// default explicit
template <class traits>
struct Members<0, traits>
{
  static constexpr auto binding_sentinel = traits::binding_sentinel;
  using At = ::pressio::rom::impl::FomObjMixin<
    typename traits::fom_system_type, binding_sentinel>;

  using Bt = ::pressio::rom::impl::FomStatesMngrMixin<
    At,
    typename traits::fom_state_type,
    typename traits::fom_state_reconstr_type,
    typename traits::fom_states_manager_type>;

  using Ct = ProjectorMixin<Bt, typename traits::projector_type>;
  using Dt = DefaultExplicitSystemMixin<Ct, typename traits::rom_system_type>;
  using type = ::pressio::rom::impl::ExplicitStepperMixin<
    Dt, typename traits::stepper_type>;
};

// default implicit
template <class traits>
struct Members<1, traits>
{
  static constexpr auto binding_sentinel = traits::binding_sentinel;
  using At = ::pressio::rom::impl::FomObjMixin<
    typename traits::fom_system_type, binding_sentinel>;

  using Bt = ::pressio::rom::impl::FomStatesMngrMixin<
    At,
    typename traits::fom_state_type,
    typename traits::fom_state_reconstr_type,
    typename traits::fom_states_manager_type>;

  using Ct  = ProjectorMixin<Bt, typename traits::projector_type>;
  using Dt  = DefaultImplicitPoliciesMixin<
    Ct,
    typename traits::residual_policy_type,
    typename traits::jacobian_policy_type>;
  using type = ::pressio::rom::impl::ImplicitStepperMixin<
    Dt, typename traits::stepper_type>;
};

// masked explicit
template <class traits>
struct Members<2, traits>
{
  static constexpr auto binding_sentinel = traits::binding_sentinel;
  using At = ::pressio::rom::impl::FomObjMixin<
    typename traits::fom_system_type, binding_sentinel>;
  using Bt = ::pressio::rom::impl::FomStatesMngrMixin<
    At,
    typename traits::fom_state_type,
    typename traits::fom_state_reconstr_type,
    typename traits::fom_states_manager_type>;

  using Ct   = MaskerMixin<Bt, typename traits::masker_type>;
  using Dt   = ProjectorMixin<Ct, typename traits::projector_type>;
  using Et   = MaskedVeloExplicitSystemMixin<Dt, typename traits::rom_system_type>;
  using type = ::pressio::rom::impl::ExplicitStepperMixin<
    Et, typename traits::stepper_type>;
};

// masked implicit
template <class traits>
struct Members<3, traits>
{
  static constexpr auto binding_sentinel = traits::binding_sentinel;
  using At = ::pressio::rom::impl::FomObjMixin<
    typename traits::fom_system_type, binding_sentinel>;
  using Bt = ::pressio::rom::impl::FomStatesMngrMixin<
    At,
    typename traits::fom_state_type,
    typename traits::fom_state_reconstr_type,
    typename traits::fom_states_manager_type>;

  using Ct  = MaskerMixin<Bt, typename traits::masker_type>;
  using Dt  = ProjectorMixin<Ct, typename traits::projector_type>;
  using Et  = MaskedVeloImplicitPoliciesMixin<
    Dt, typename traits::residual_policy_type, typename traits::jacobian_policy_type>;
  using type= ::pressio::rom::impl::ImplicitStepperMixin<
    Et, typename traits::stepper_type>;
};


// hypred explicit
template <class traits>
struct Members<4, traits>
{
  static constexpr auto binding_sentinel = traits::binding_sentinel;
  using At = ::pressio::rom::impl::FomObjMixin<
    typename traits::fom_system_type, binding_sentinel>;
  using Bt = ::pressio::rom::impl::FomStatesMngrMixin<
    At,
    typename traits::fom_state_type,
    typename traits::fom_state_reconstr_type,
    typename traits::fom_states_manager_type>;

  using Ct   = ProjectorMixin<Bt, typename traits::projector_type>;
  using Dt = HypRedVeloExplicitSystemMixin<Ct, typename traits::rom_system_type>;
  using type = ::pressio::rom::impl::ExplicitStepperMixin<
    Dt, typename traits::stepper_type>;
};

// hypred implicit
template <class traits>
struct Members<5, traits>
{
  static constexpr auto binding_sentinel = traits::binding_sentinel;
  using At = ::pressio::rom::impl::FomObjMixin<
    typename traits::fom_system_type, binding_sentinel>;
  using Bt = ::pressio::rom::impl::FomStatesMngrMixin<
    At,
    typename traits::fom_state_type,
    typename traits::fom_state_reconstr_type,
    typename traits::fom_states_manager_type>;

  using Ct  = ProjectorMixin<Bt, typename traits::projector_type>;
  using Dt  = HypRedVeloImplicitPoliciesMixin<
    Ct, typename traits::residual_policy_type, typename traits::jacobian_policy_type>;
  using type= ::pressio::rom::impl::ImplicitStepperMixin<
    Dt, typename traits::stepper_type>;
};

template <int flag, typename ...Args>
class ProblemContinuousTimeApi
{
public:
  using traits = ::pressio::Traits<ProblemContinuousTimeApi<flag, Args...>>;

private:
  using fom_system_type	     = typename traits::fom_system_type;
  using decoder_type	     = typename traits::decoder_type;
  using galerkin_state_type     = typename traits::galerkin_state_type;
  using stepper_type	     = typename traits::stepper_type;
  using fom_state_type	     = typename traits::fom_state_type;
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
  ProblemContinuousTimeApi() = delete;
  ProblemContinuousTimeApi(const ProblemContinuousTimeApi &) = default;
  ProblemContinuousTimeApi & operator=(const ProblemContinuousTimeApi &) = delete;
  ProblemContinuousTimeApi(ProblemContinuousTimeApi &&) = default;
  ProblemContinuousTimeApi & operator=(ProblemContinuousTimeApi &&) = delete;
  ~ProblemContinuousTimeApi() = default;

  template<
    int _flag = flag,
    mpl::enable_if_t<_flag==0 or _flag==1, int> = 0
    >
  ProblemContinuousTimeApi(const fom_system_type & fomObj,
			   decoder_type & decoder,
			   const galerkin_state_type & romState,
			   const fom_state_type & fomNominalState)
    : members_(romState, fomObj, decoder, fomNominalState){}

  template<
    int _flag = flag, class ...Args2,
    mpl::enable_if_t<_flag >=2 and _flag <= 5, int> = 0
    >
  ProblemContinuousTimeApi(const fom_system_type & fomObj,
			   decoder_type & decoder,
			   const galerkin_state_type & romState,
			   const fom_state_type & fomNominalState,
			   const typename traits::projector_type & projector,
			   Args2 && ...args)
    : members_(romState, fomObj, decoder, fomNominalState,
	       projector, std::forward<Args2>(args) ...){}
};

}}}}//end namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_ROM_GALERKIN_DEFAULT_PROBLEM_EXPLICIT_STEPPING_HPP_
