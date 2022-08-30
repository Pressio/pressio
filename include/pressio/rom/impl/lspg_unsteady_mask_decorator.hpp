/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_cont_time_decorators.hpp
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

#ifndef ROM_IMPL_ROM_LSPG_UNSTEADY_CONT_TIME_MASK_DECORATOR_HPP_
#define ROM_IMPL_ROM_LSPG_UNSTEADY_CONT_TIME_MASK_DECORATOR_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  class ResidualMaskerType,
  class JacobianMaskerType,
  class Maskable
  >
class LspgMaskDecorator : public Maskable
{
  using unmasked_residual_type = typename Maskable::residual_type;
  using unmasked_jacobian_type = typename Maskable::jacobian_type;
  using ind_var_type = typename Maskable::independent_variable_type;

public:
  // required
  using independent_variable_type = typename Maskable::independent_variable_type;
  using state_type    = typename Maskable::state_type;
  using residual_type = typename ResidualMaskerType::result_type;
  using jacobian_type = typename JacobianMaskerType::result_type;

  LspgMaskDecorator() = delete;

  template <
    class TrialSpaceType,
    class FomSystemType,
    template<class> class LspgFomStatesManager
    >
  LspgMaskDecorator(const TrialSpaceType & trialSpace,
		    const FomSystemType & fomSystem,
		    LspgFomStatesManager<TrialSpaceType> & fomStatesManager,
		    const ResidualMaskerType & rMasker,
		    const JacobianMaskerType & jaMasker)
    : Maskable(trialSpace, fomSystem, fomStatesManager),
      rMasker_(rMasker),
      jMasker_(jaMasker),
      unMaskedResidual_(Maskable::createResidual()),
      unMaskedJacobian_(Maskable::createJacobian())
  {}

public:
  state_type createState() const{
    return Maskable::createState();
  }

  residual_type createResidual() const{
    auto tmp = Maskable::createResidual();
    ::pressio::ops::set_zero(tmp);
    return rMasker_.get().createApplyMaskResult(tmp);
  }

  jacobian_type createJacobian() const{
    auto tmp = Maskable::createJacobian();
    ::pressio::ops::set_zero(tmp);
    return jMasker_.get().createApplyMaskResult(tmp);
  }

  template <class StencilStatesContainerType, class StencilRhsContainerType>
  void operator()(::pressio::ode::StepScheme odeSchemeName,
		  const state_type & predictedReducedState,
		  const StencilStatesContainerType & reducedStatesStencilManager,
		  StencilRhsContainerType & fomRhsStencilManger,
		  const ::pressio::ode::StepEndAt<ind_var_type> & rhsEvaluationTime,
		  ::pressio::ode::StepCount step,
		  const ::pressio::ode::StepSize<ind_var_type> & dt,
		  residual_type & R,
		  jacobian_type & J,
		  bool computeJacobian) const
  {
    Maskable::operator()(odeSchemeName, predictedReducedState,
			 reducedStatesStencilManager, fomRhsStencilManger,
			 rhsEvaluationTime, step, dt,
			 unMaskedResidual_, unMaskedJacobian_,
			 computeJacobian);
    rMasker_(unMaskedResidual_, R);
    jMasker_(unMaskedJacobian_, J);
  }

private:
  mutable unmasked_residual_type unMaskedResidual_;
  mutable unmasked_jacobian_type unMaskedJacobian_;
  std::reference_wrapper<const ResidualMaskerType> rMasker_;
  std::reference_wrapper<const JacobianMaskerType> jMasker_;
};

}}}
#endif  // ROM_IMPL_ROM_LSPG_UNSTEADY_CONT_TIME_DECORATORS_HPP_