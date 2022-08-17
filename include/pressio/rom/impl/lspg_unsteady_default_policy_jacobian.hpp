/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_policy_residual.hpp
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

#ifndef ROM_IMPL_ROM_LSPG_UNSTEADY_DEFAULT_POLICY_JACOBIAN_HPP_
#define ROM_IMPL_ROM_LSPG_UNSTEADY_DEFAULT_POLICY_JACOBIAN_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  class IndVarType,
  class ReducedStateType,
  class LspgJacobianType,
  class TrialSpaceType,
  class FomSystemType
  >
class LspgUnsteadyJacobianPolicy
{
public:
  // required
  using independent_variable_type = IndVarType;
  using state_type    = ReducedStateType;
  using jacobian_type = LspgJacobianType;

public:
  LspgUnsteadyJacobianPolicy() = delete;
  LspgUnsteadyJacobianPolicy(const TrialSpaceType & trialSpace,
			     const FomSystemType & fomSystem,
			     LspgFomStatesManager<TrialSpaceType> & fomStatesManager)
    : trialSpace_(trialSpace),
      fomSystem_(fomSystem),
      fomStatesManager_(fomStatesManager)
  {}

public:
  state_type createState() const{
    // this needs to create an instance of the reduced state
    return trialSpace_.get().createReducedState();
  }

  jacobian_type createJacobian() const{
    // lspg jacobian is of the same shape as the basis
    auto J = ::pressio::ops::clone(trialSpace_.get().viewBasis());
    ::pressio::ops::set_zero(J);
    return J;
  }

  template <class StencilStatesContainerType>
  void operator()(::pressio::ode::StepScheme odeSchemeName,
		  const state_type & predictedReducedState,
		  const StencilStatesContainerType & /*unused*/,
		  const ::pressio::ode::StepEndAt<IndVarType> & fomEvaluationTime,
		  ::pressio::ode::StepCount /*step*/,
		  const ::pressio::ode::StepSize<IndVarType> & dt,
		  jacobian_type & J) const
  {

    const auto & phi = trialSpace_.get().viewBasis();
    const auto & fomState = fomStatesManager_(::pressio::ode::nPlusOne());
    fomSystem_.get().applyJacobian(fomState, phi, fomEvaluationTime.get(), J);

    // for lspg, e.g. BDF1,2, we have something like:
    //
    //	 lspgJac = decoderJac + dt*coeff*J*decoderJac
    //
    // where J is the fom jacobian and coeff depends on the scheme.
    //
    // Here, J already contains J*phi (see above),
    // so we just need to update J properly.

    using basis_sc_t = typename ::pressio::Traits<
      typename TrialSpaceType::basis_type>::scalar_type;
    const auto one = ::pressio::utils::Constants<basis_sc_t>::one();

    if (odeSchemeName == ::pressio::ode::StepScheme::BDF1){
      const auto factor = dt.get()*::pressio::ode::constants::bdf1<IndVarType>::c_f_;
      ::pressio::ops::update(J, factor, phi, one);
    }
    else{
      throw std::runtime_error("Only BDF1 currently impl for default unstedy LSPG");
    }
  }

private:
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  std::reference_wrapper<LspgFomStatesManager<TrialSpaceType>> fomStatesManager_;
};

}}}
#endif  // ROM_IMPL_ROM_LSPG_UNSTEADY_POLICY_JACOBIAN_HPP_
