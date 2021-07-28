/*
//@HEADER
// ************************************************************************
//
// solvers_hessian_gradient_corrector.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_CORRECTION_MIXINS_SOLVERS_HESSIAN_GRADIENT_CORRECTOR_HPP_
#define SOLVERS_NONLINEAR_IMPL_CORRECTION_MIXINS_SOLVERS_HESSIAN_GRADIENT_CORRECTOR_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template<typename T, typename state_type, typename lin_solver_t>
class HessianGradientCorrector : public T
{
public:
  using state_t = state_type;
  using sc_t = typename state_type::traits::scalar_t;
  using state_wrapped_t = typename state_type::traits::wrapped_t;

private:
  state_type correction_ = {};
  ::pressio::utils::instance_or_reference_wrapper<lin_solver_t> solverObj_;
  sc_t residNormCurrCorrStep_ = {};
  sc_t gradientNormCurrCorrStep_ = {};
  sc_t correctionNormCurrCorrStep_ = {};

public:
  HessianGradientCorrector() = delete;

  template <typename system_t, typename lsT, typename ...Args>
  HessianGradientCorrector(const system_t & system,
			   const state_type & state,
			   lsT && solverIn,
			   Args && ... args)
    : T(system, state, std::forward<Args>(args)...),
      correction_(state),
      solverObj_(std::forward<lsT>(solverIn))
  {
    constexpr auto zero = ::pressio::utils::constants<sc_t>::zero();
    ::pressio::ops::fill(correction_, zero);
  }

  template <typename system_t, typename lsT, typename ...Args>
  HessianGradientCorrector(const system_t & system,
			   const state_wrapped_t & state,
			   lsT && solverIn,
			   Args && ... args)
    : HessianGradientCorrector(system,
			       state_type(state) /*needs a wrapped object*/,
			       std::forward<lsT>(solverIn),
			       std::forward<Args>(args)...)
  {}

  HessianGradientCorrector(HessianGradientCorrector const &) = default;
  HessianGradientCorrector & operator=(HessianGradientCorrector const &) = default;
  HessianGradientCorrector(HessianGradientCorrector && o) = default;
  HessianGradientCorrector & operator=(HessianGradientCorrector && o) = default;
  ~HessianGradientCorrector() = default;

public:
  template <typename system_t>
  void computeCorrection(const system_t & sys,
			 state_type & state,
			 bool recomputeSystemJacobian = true)
  {
    PRESSIOLOG_DEBUG("hessian/gradient correction");
    T::computeOperators(sys, state,
			residNormCurrCorrStep_,
			recomputeSystemJacobian);

    const auto & H = T::hessianCRef();
    const auto & g = T::gradientCRef();
    solverObj_.get().solve(H, g, correction_);

    correctionNormCurrCorrStep_ = pressio::ops::norm2(correction_);
    gradientNormCurrCorrStep_ = pressio::ops::norm2(g);
  }

  bool hasGradientComputation() const {
    return true;
  }

  void resetForNewCall(){
    T::resetForNewCall();
  }

  const state_type & correctionCRef() const{
    return correction_;
  }

  const sc_t & correctionNormCurrentCorrectionStep() const{
    return correctionNormCurrCorrStep_;
  }

  const sc_t & gradientNormCurrentCorrectionStep() const{
    return gradientNormCurrCorrStep_;
  }

  const sc_t & residualNormCurrentCorrectionStep() const{
    return residNormCurrCorrStep_;
  }
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_CORRECTION_MIXINS_SOLVERS_HESSIAN_GRADIENT_CORRECTOR_HPP_
