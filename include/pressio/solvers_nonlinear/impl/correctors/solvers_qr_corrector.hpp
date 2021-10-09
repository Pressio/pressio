/*
//@HEADER
// ************************************************************************
//
// solvers_qr_corrector.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_CORRECTORS_SOLVERS_QR_CORRECTOR_HPP_
#define SOLVERS_NONLINEAR_IMPL_CORRECTORS_SOLVERS_QR_CORRECTOR_HPP_

namespace pressio{ namespace nonlinearsolvers{ namespace impl{

template<class T, class StateType, class QRSolverType>
class QRCorrector : public T
{
public:
  using typename T::scalar_type;
  using state_type = StateType;

private:
  state_type correction_ = {};
  state_type QTResid_ = {};
  state_type g_ = {};

  ::pressio::utils::InstanceOrReferenceWrapper<QRSolverType> solverObj_;
  scalar_type residNormCurrCorrStep_ = {};
  scalar_type gradientNormCurrCorrStep_ = {};
  scalar_type correctionNormCurrCorrStep_ = {};

public:
  QRCorrector() = delete;

  template <typename SystemType, typename qrsT>
  QRCorrector(const SystemType & system,
	      const state_type & state,
	      qrsT && solverObj)
    : T(system, state),
      correction_(::pressio::ops::clone(state)),
      QTResid_(::pressio::ops::clone(state)),
      g_(::pressio::ops::clone(state)),
      solverObj_(std::forward<qrsT>(solverObj))
  {
    constexpr auto zero = ::pressio::utils::Constants<scalar_type>::zero();
    ::pressio::ops::fill(correction_, zero);
    ::pressio::ops::fill(QTResid_, zero);
    ::pressio::ops::fill(g_, zero);
  }

  QRCorrector(QRCorrector const &) = default;
  QRCorrector & operator=(QRCorrector const &) = default;
  QRCorrector(QRCorrector && o) = default;
  QRCorrector & operator=(QRCorrector && o) = default;
  ~QRCorrector() = default;

public:
  template <typename SystemType>
  void computeCorrection(const SystemType & sys,
			 state_type & state,
			 bool recomputeSystemJacobian = true)
  {
    PRESSIOLOG_DEBUG("QR-based correction");
    T::computeOperators(sys, state,
			residNormCurrCorrStep_,
			recomputeSystemJacobian);
    const auto & r = T::residualCRef();
    const auto & J = T::jacobianCRef();

    // J = QR
    solverObj_.get().computeThin(J);
    // compute Q^T r
    solverObj_.get().applyQTranspose(r, QTResid_);
    // compute gradient = R^T Q^T r
    solverObj_.get().applyRTranspose(QTResid_, g_);

    // solve: R correction = Q^T Residual
    solverObj_.get().solve(QTResid_, correction_);
    // scale by -1 for sign convention
    pressio::ops::scale(correction_, utils::Constants<scalar_type>::negOne() );

    correctionNormCurrCorrStep_ = pressio::ops::norm2(correction_);
    gradientNormCurrCorrStep_ = pressio::ops::norm2(g_);
  }

  bool hasGradientComputation() const { return true; }

  void resetForNewCall(){
    T::resetForNewCall();
  }

  const state_type & correctionCRef() const{
    return correction_;
  }

  const state_type & gradientCRef() const{
    return g_;
  }

  const scalar_type & correctionNormCurrentCorrectionStep() const{
    return correctionNormCurrCorrStep_;
  }

  const scalar_type & gradientNormCurrentCorrectionStep() const{
    return gradientNormCurrCorrStep_;
  }

  const scalar_type & residualNormCurrentCorrectionStep() const{
    return residNormCurrCorrStep_;
  }

};

}}}
#endif  // SOLVERS_NONLINEAR_IMPL_CORRECTORS_SOLVERS_QR_CORRECTOR_HPP_
