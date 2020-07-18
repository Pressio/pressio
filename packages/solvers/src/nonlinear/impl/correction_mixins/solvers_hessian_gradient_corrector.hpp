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

template<
  typename T, typename state_t, typename lin_solver_t, ::pressio::Norm normType
  >
class HessianGradientCorrector : public T
{
  using sc_t = typename ::pressio::containers::details::traits<state_t>::scalar_t;

  state_t correction_ = {};

  typename std::conditional<
    ::pressio::ops::predicates::is_object_pybind<lin_solver_t>::value,
    lin_solver_t, lin_solver_t &>::type solverObj_;

  sc_t residualNorm_ = {};
  sc_t gradientNorm_ = {};
  sc_t correctionNorm_ = {};

public:
  static constexpr auto normType_ = normType;

  HessianGradientCorrector() = delete;

  template <typename system_t, typename ...Args>
  HessianGradientCorrector(const system_t & system,
			   const state_t & state,
			   lin_solver_t & solverObj,
			   Args && ... args)
    : T(system, state, std::forward<Args>(args)...),
      correction_(state),
      solverObj_(solverObj)
  {}

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  // this is the overload use for pressio4py since the state object passed to the constructor
  // is not a wrapper bbut a native pybind array, so we need to use the wrapped type.
  // However, this is only true for the constructor. For the actual computation below,
  // the state is wrapper type.
  template <typename system_t, typename ...Args>
  HessianGradientCorrector(const system_t & system,
			   const typename ::pressio::containers::details::traits<state_t>::wrapped_t & state,
			   lin_solver_t & solverObj,
			   Args && ... args)
    : T(system, state_t(state) /*the parent needs a wrapped object*/, std::forward<Args>(args)...),
      correction_(state),
      solverObj_(solverObj)
  {}
#endif

public:
  template <typename system_t>
  void computeCorrection(const system_t & sys, state_t & state)
  {
    T::computeOperators(sys, state, normType, residualNorm_);
    auto & H = T::getHessian();
    auto & g = T::getGradient();
    this->doLinearSolve(H, g);

    correctionNorm_ = pressio::ops::norm2(correction_);
    gradientNorm_ = pressio::ops::norm2(g);
  }

  const state_t & viewCorrection() const{ return correction_; }
  const sc_t correctionNormCurrentCorrectionStep() const{ return correctionNorm_; }
  const sc_t gradientNormCurrentCorrectionStep() const{ return gradientNorm_; }
  const sc_t residualNormCurrentCorrectionStep() const{ return residualNorm_; }

  template< typename system_t>
  void residualNorm(const system_t & system, const state_t & state, sc_t & result){
    T::residualNorm(system, state, normType, result);
  }

private:
  template <typename H_t, typename g_t, typename _lin_solver_t = lin_solver_t>
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  mpl::enable_if_t<!std::is_same<_lin_solver_t, pybind11::object>::value>
#else
  void
#endif
  doLinearSolve(H_t & H, const g_t & g)
  {
    solverObj_.solveAllowMatOverwrite(H, g, correction_);
  }

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  template <typename H_t, typename g_t, typename _lin_solver_t = lin_solver_t>
  mpl::enable_if_t<std::is_same<_lin_solver_t, pybind11::object>::value>
  doLinearSolve(H_t & H, const g_t & g)
  {
    solverObj_.attr("solve")(*H.data(), *g.data(), *correction_.data());
  }
#endif

};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_CORRECTION_MIXINS_SOLVERS_HESSIAN_GRADIENT_CORRECTOR_HPP_
