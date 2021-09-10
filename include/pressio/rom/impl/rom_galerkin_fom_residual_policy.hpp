/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_fom_residual_policy.hpp
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

#ifndef ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_FOM_RESIDUAL_POLICY_HPP_
#define ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_FOM_RESIDUAL_POLICY_HPP_

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{

template <class FomStatesManagerType, class FomResidualType, class FomSystemType>
class FomResidualPolicyDiscreteTimeApi
{
public:
  FomResidualPolicyDiscreteTimeApi() = delete;
  FomResidualPolicyDiscreteTimeApi(const FomResidualPolicyDiscreteTimeApi &) = default;
  FomResidualPolicyDiscreteTimeApi & operator=(const FomResidualPolicyDiscreteTimeApi &) = delete;
  FomResidualPolicyDiscreteTimeApi(FomResidualPolicyDiscreteTimeApi &&) = default;
  FomResidualPolicyDiscreteTimeApi & operator=(FomResidualPolicyDiscreteTimeApi &&) = delete;
  ~FomResidualPolicyDiscreteTimeApi() = default;

  FomResidualPolicyDiscreteTimeApi(const FomSystemType & fomSystem,
				   FomStatesManagerType & fomStatesMngr)
    : fomStatesMngr_(fomStatesMngr),
      fomSystem_(fomSystem),
      fomResidual_(fomSystem.createDiscreteTimeResidual())
  {
    ::pressio::ops::set_zero(fomResidual_);
  }

public:
  const FomResidualType & get() const{
    return fomResidual_;
  }

  template <class GalerkinStateType, class ScalarType>
  void compute(const GalerkinStateType & galerkin_state_np1,
	       const ScalarType & time_np1,
	       const ScalarType & dt,
	       const int32_t & step_number) const
  {
    this->doFomStatesReconstruction(step_number, galerkin_state_np1);

    const auto & fom_state_np1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    fomSystem_.get().discreteTimeResidual(step_number, time_np1, dt,
				    fomResidual_, fom_state_np1);
  }

  template <class GalerkinStateType, class ScalarType>
  void compute(const GalerkinStateType & galerkin_state_np1,
	       const ScalarType & time_np1,
	       const ScalarType & dt,
	       const int32_t & step_number,
	       const GalerkinStateType & galerkin_state_n) const
  {
    this->doFomStatesReconstruction(step_number, galerkin_state_np1, galerkin_state_n);

    const auto & fom_state_np1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    const auto & fom_state_n   = fomStatesMngr_(::pressio::ode::n());
    fomSystem_.get().discreteTimeResidual(step_number, time_np1, dt,
				    fomResidual_, fom_state_np1, fom_state_n);
  }

  template <class GalerkinStateType, class ScalarType>
  void compute(const GalerkinStateType & galerkin_state_np1,
	       const ScalarType & time_np1,
	       const ScalarType & dt,
	       const int32_t & step_number,
	       const GalerkinStateType & galerkin_state_n,
	       const GalerkinStateType & galerkin_state_nm1) const
  {
    this->doFomStatesReconstruction(step_number, galerkin_state_np1,
				    galerkin_state_n, galerkin_state_nm1);

    const auto & fom_state_np1 = fomStatesMngr_(::pressio::ode::nPlusOne());
    const auto & fom_state_n   = fomStatesMngr_(::pressio::ode::n());
    const auto & fom_state_nm1 = fomStatesMngr_(::pressio::ode::nMinusOne());

    fomSystem_.get().discreteTimeResidual(step_number, time_np1, dt, fomResidual_,
				      fom_state_np1, fom_state_n, fom_state_nm1);
  }

private:
  template <class GalerkinStateType>
  void doFomStatesReconstruction(const int32_t & step_number,
				 const GalerkinStateType & galerkin_state_np1) const
  {
    fomStatesMngr_.get().reconstructAt(galerkin_state_np1, ::pressio::ode::nPlusOne());
  }

  template <class GalerkinStateType>
  void doFomStatesReconstruction(const int32_t & step_number,
				 const GalerkinStateType & galerkin_state_np1,
				 const GalerkinStateType & galerkin_state_n) const
  {
    /* the FOM state corresponding to the new predicted state has to be
     * recomputed every time regardless of the time step chaning or not,
     *  since we might be inside a non-linear solve
     * where the time step does not change but this residual method
     * is called multiple times. */
    fomStatesMngr_.get().reconstructAt(galerkin_state_np1, ::pressio::ode::nPlusOne());

    /* previous FOM states should only be recomputed when the time step changes.
     * The method below does not recompute all previous states, but only
     * recomputes the n-th state and updates/shifts back all the other
     * FOM states stored. */
    if (storedStep_ != step_number){
      fomStatesMngr_.get().reconstructAtAndUpdatePrevious(galerkin_state_n,
							  ::pressio::ode::n());
      storedStep_ = step_number;
    }
  }

  template <class GalerkinStateType>
  void doFomStatesReconstruction(const int32_t & step_number,
				 const GalerkinStateType & galerkin_state_np1,
				 const GalerkinStateType & galerkin_state_n,
				 const GalerkinStateType & galerkin_state_nm1) const
  {
    (void)galerkin_state_nm1;
    doFomStatesReconstruction(step_number, galerkin_state_np1, galerkin_state_n);
  }

private:
  // storedStep is used to keep track of which step we are doing.
  // This is used to decide whether we need to update/recompute the previous
  // FOM states or not. Since it does not make sense to recompute previous
  // FOM states if we are not in a new time step.
  mutable int32_t storedStep_ = {};

  std::reference_wrapper<FomStatesManagerType> fomStatesMngr_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable FomResidualType fomResidual_;
};

}}}}//end namespace
#endif  // ROM_GALERKIN_IMPL_DISCRETE_TIME_API_POLICIES_ROM_GALERKIN_FOM_RESIDUAL_POLICY_HPP_
