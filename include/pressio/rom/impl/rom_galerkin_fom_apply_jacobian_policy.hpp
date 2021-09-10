/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_fom_apply_jacobian_policy.hpp
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

#ifndef ROM_GALERKIN_IMPL_FOM_APPLY_JACOBIAN_POLICY_HPP_
#define ROM_GALERKIN_IMPL_FOM_APPLY_JACOBIAN_POLICY_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <
  bool is_continuous_time,
  class FomStatesManagerType,
  class fom_apply_jac_type,
  class DecoderType,
  class FomSystemType
  >
class DefaultFomApplyJacobianEvaluator
{

public:
  DefaultFomApplyJacobianEvaluator() = delete;
  DefaultFomApplyJacobianEvaluator(const DefaultFomApplyJacobianEvaluator &) = default;
  DefaultFomApplyJacobianEvaluator & operator=(const DefaultFomApplyJacobianEvaluator &) = delete;
  DefaultFomApplyJacobianEvaluator(DefaultFomApplyJacobianEvaluator &&) = default;
  DefaultFomApplyJacobianEvaluator & operator=(DefaultFomApplyJacobianEvaluator &&) = delete;
  ~DefaultFomApplyJacobianEvaluator() = default;

  template<
    bool _is_cont_time = is_continuous_time,
    mpl::enable_if_t<_is_cont_time, int> = 0
    >
  DefaultFomApplyJacobianEvaluator(const FomSystemType & fomSystem,
				   FomStatesManagerType & fomStatesMngr,
				   const DecoderType & decoder)
    : fomStatesMngr_(fomStatesMngr),
      fomSystem_(fomSystem),
      phi_(decoder.jacobianCRef()),
      fomApplyJac_(fomSystem.createApplyJacobianResult(phi_.get()))
  {
    ::pressio::ops::set_zero(fomApplyJac_);
  }

  template<
    bool _is_cont_time = is_continuous_time,
    mpl::enable_if_t<!_is_cont_time, int> = 0
    >
  DefaultFomApplyJacobianEvaluator(const FomSystemType & fomSystem,
				   FomStatesManagerType & fomStatesMngr,
				   const DecoderType & decoder)
    : fomStatesMngr_(fomStatesMngr),
      fomSystem_(fomSystem),
      phi_(decoder.jacobianCRef()),
      fomApplyJac_(fomSystem.createApplyDiscreteTimeJacobianResult(phi_.get()))
  {
    ::pressio::ops::set_zero(fomApplyJac_);
  }

public:
  const fom_apply_jac_type & get() const{ return fomApplyJac_; }

  template<class GalerkinStateType, class ScalarType, bool _is_cont_time = is_continuous_time>
  mpl::enable_if_t<_is_cont_time>
  compute(const GalerkinStateType & galerkinState,
	  const ScalarType & evaluationTime) const
  {
    /* here we assume that current FOM state has been reconstructd
       by the residual policy. So we do not recompute the FOM state.
       Maybe we should find a way to ensure this is the case.
     */
    const auto & fomState = fomStatesMngr_(::pressio::ode::nPlusOne());

    // call applyJacobian on the fom object
    fomSystem_.get().applyJacobian(fomState, phi_.get(), evaluationTime, fomApplyJac_);
  }


  template<class galerkin_state_t,  class scalar_t, bool _is_cont_time = is_continuous_time>
  mpl::enable_if_t< !_is_cont_time >
  compute(const galerkin_state_t & galerkin_state_np1,
	  const scalar_t & timeAtNextStep,
	  const scalar_t & dt,
	  const int32_t & currentStepNumber) const
  {
    // we don't use gal_states because we assume the residual policy
    // already used them to reconstuct the FOM states
    // (to fix that at some point)

    const auto & ynp1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    fomSystem_.get().applyDiscreteTimeJacobian(currentStepNumber, timeAtNextStep,
					   dt, phi_.get(),
					   fomApplyJac_,
					   ynp1);
  }

  template<class galerkin_state_t,  class scalar_t, bool _is_cont_time = is_continuous_time>
  mpl::enable_if_t< !_is_cont_time >
  compute(const galerkin_state_t & galerkin_state_np1,
	  const scalar_t & timeAtNextStep,
	  const scalar_t & dt,
	  const int32_t & currentStepNumber,
	  const galerkin_state_t & galerkin_state_n) const
  {
    // we don't use gal_states because we assume the residual policy
    // already used them to reconstuct the FOM states
    // (to fix that at some point)

    const auto & ynp1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    const auto & yn   = fomStatesMngr_(::pressio::ode::n());
    fomSystem_.get().applyDiscreteTimeJacobian(currentStepNumber, timeAtNextStep,
					   dt, phi_.get(),
					   fomApplyJac_,
					   ynp1, yn);
  }

  template<class galerkin_state_t,  class scalar_t, bool _is_cont_time = is_continuous_time>
  mpl::enable_if_t< !_is_cont_time >
  compute(const galerkin_state_t & galerkin_state_np1,
	  const scalar_t & timeAtNextStep,
	  const scalar_t & dt,
	  const int32_t & currentStepNumber,
	  const galerkin_state_t & galerkin_state_n,
	  const galerkin_state_t & galerkin_state_nm1) const
  {
    // we don't use gal_states because we assume the residual policy
    // already used them to reconstuct the FOM states
    // (to fix that at some point)

    const auto & ynp1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    const auto & yn   = fomStatesMngr_(::pressio::ode::n());
    const auto & ynm1 = fomStatesMngr_(::pressio::ode::nMinusOne());
    fomSystem_.get().applyDiscreteTimeJacobian(currentStepNumber, timeAtNextStep,
					   dt, phi_.get(),
					   fomApplyJac_,
					   ynp1, yn);
  }

private:
  std::reference_wrapper<FomStatesManagerType> fomStatesMngr_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  std::reference_wrapper<const typename DecoderType::jacobian_type> phi_;
  mutable fom_apply_jac_type fomApplyJac_ = {};

};

}}}}//end namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_CONTINUOUS_TIME_API_POLICIES_ROM_GALERKIN_FOM_APPLY_JACOBIAN_POLICY_HPP_
