/*
//@HEADER
// ************************************************************************
//
// solvers_hessian_gradient_corrector_impl.hpp
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

#ifndef PRESSIO_SOLVERS_HESSIAN_GRADIENT_CORRECTOR_IMPL_HPP_
#define PRESSIO_SOLVERS_HESSIAN_GRADIENT_CORRECTOR_IMPL_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template<
  typename T, typename state_t, typename lin_solver_t, ::pressio::solvers::Norm normType
  >
class HessianGradientCorrector : public T
{
  using sc_t = typename ::pressio::containers::details::traits<state_t>::scalar_t;

  state_t correction_ = {};
  lin_solver_t & solverObj_;
  sc_t residualNorm_ = {};

public:
  static constexpr auto normType_ = normType;

  HessianGradientCorrector() = delete;

  template <typename system_t, typename ...Args>
  HessianGradientCorrector(const system_t & system, const state_t & state, lin_solver_t & solverObj)
    : T(system, state), correction_(state), solverObj_(solverObj){}

public:
  template <typename system_t>
  void computeCorrection(const system_t & sys, state_t & state)
  {
    T::computeOperators(sys, state, normType, residualNorm_);
    auto & H = T::getHessian();
    auto & g = T::getGradient();
    solverObj_.solve(H, g, correction_);
    std::cout << std::fixed
    	      << std::setprecision(15)
    	      << residualNorm_ << " "
    	      << pressio::ops::norm2(g) << " "
    	      << pressio::ops::norm2(correction_)
    	      << std::endl;
  }

  const state_t & viewCorrection() const{ return correction_; }
  const sc_t currentResidualNorm() const{ return residualNorm_; }

  template< typename system_t>
  void residualNorm(const system_t & system, const state_t & state, sc_t & result) const{
    system.residualNorm(state, normType, result);
  }
};

}}}}
#endif
