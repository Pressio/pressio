/*
//@HEADER
// ************************************************************************
//
// solvers_create_public_api.hpp
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

#ifndef SOLVERS_NONLINEAR_SOLVERS_CREATE_GAUSS_NEWTON_HPP_
#define SOLVERS_NONLINEAR_SOLVERS_CREATE_GAUSS_NEWTON_HPP_

#include "./impl/solvers_nonlinear_compose.hpp"

namespace pressio{ namespace nonlinearsolvers{

// ----------------------------------------------------------------
// Gauss-Newton
// ----------------------------------------------------------------
template<class SystemType, class LinearSolverType>
#ifdef PRESSIO_ENABLE_CXX20
requires
     (OverdeterminedSystemWithResidualAndJacobian<SystemType>
   || OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>
   || SystemWithHessianAndGradient<SystemType>
   || SystemWithFusedHessianAndGradient<SystemType>)
  && LinearSolverForNonlinearLeastSquares<
       mpl::remove_cvref_t<LinearSolverType>,
       typename SystemType::state_type>
#endif
auto create_gauss_newton(const SystemType & system,
			 LinearSolverType && linSolver)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert
    (   OverdeterminedSystemWithResidualAndJacobian<SystemType>::value
     or OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>::value
     or SystemWithHessianAndGradient<SystemType>::value
     or SystemWithFusedHessianAndGradient<SystemType>::value,
     "Gauss-Newton: system does not satisfy the residual/jacobian, or hessian/gradient concepts.");

  static_assert(LinearSolverForNonlinearLeastSquares<
		mpl::remove_cvref_t<LinearSolverType>,
		typename SystemType::state_type>::value,
		"Gauss-newton: linear solver not satisfying the concept.");
#endif

  return impl::ComposeGaussNewton_t<SystemType, LinearSolverType>
    (system, std::forward<LinearSolverType>(linSolver));
}

// ----------------------------------------------------------------
// Weighted Gauss-Newton
// ----------------------------------------------------------------

template<class SystemType, class LinearSolverType, class WeightingOpType>
#ifdef PRESSIO_ENABLE_CXX20
requires
  (OverdeterminedSystemWithResidualAndJacobian<SystemType>
   || OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>)
  && LinearSolverForNonlinearLeastSquares<
       mpl::remove_cvref_t<LinearSolverType>,
       typename SystemType::state_type>
  && LeastSquaresWeightingOperator<
       mpl::remove_cvref_t<WeightingOpType>,
       typename SystemType::residual_type,
       typename SystemType::jacobian_type>
#endif
auto create_gauss_newton(const SystemType & system,
			 LinearSolverType && linSolver,
                         WeightingOpType && weightOperator)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert
    (   OverdeterminedSystemWithResidualAndJacobian<SystemType>::value
     or OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>::value,
     "Weighted Gauss-Newton: system not satisfying the residual/jacobian concept.");

  static_assert(LinearSolverForNonlinearLeastSquares<
		mpl::remove_cvref_t<LinearSolverType>,
		typename SystemType::state_type>::value,
		"Weighted Gauss-newton: linear solver not satisfying the concept.");
  static_assert(LeastSquaresWeightingOperator<
		mpl::remove_cvref_t<WeightingOpType>,
		typename SystemType::residual_type,
		typename SystemType::jacobian_type>::value,
		"Invalid weighting operator for GN or LM");
#endif

  return impl::ComposeGaussNewton_t<SystemType, LinearSolverType, WeightingOpType>
    (system, std::forward<LinearSolverType>(linSolver),
     std::forward<WeightingOpType>(weightOperator));
}

// ----------------------------------------------------------------
// Gauss-Newton QR
// ----------------------------------------------------------------

template<class SystemType, class QRSolverType>
#ifdef PRESSIO_ENABLE_CXX20
requires
    (OverdeterminedSystemWithResidualAndJacobian<SystemType>
  || OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>)
  && QRSolverForGnQr<
	mpl::remove_cvref_t<QRSolverType>,
	typename SystemType::state_type,
	typename SystemType::jacobian_type,
	typename SystemType::residual_type
     >
#endif
auto create_gauss_newtonQR(const SystemType & system,
			   QRSolverType && qrSolver)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert
    (   OverdeterminedSystemWithResidualAndJacobian<SystemType>::value
     or OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>::value,
     "Gauss-Newton with QR: system not satisfying the residual/jacobian concept.");

  static_assert(QRSolverForGnQr<
		mpl::remove_cvref_t<QRSolverType>,
		typename SystemType::state_type,
		typename SystemType::jacobian_type,
		typename SystemType::residual_type>::value,
		"Gauss-newton QR: linear solver not satisfying the concept.");
#endif

  using state_t = typename SystemType::state_type;
  using r_t = typename SystemType::residual_type;
  using j_t = typename SystemType::jacobian_type;

  using operators_t = impl::ResidualJacobianOperators<r_t, j_t>;
  using corrector_t = impl::QRCorrector<operators_t, state_t, QRSolverType>;
  return impl::Solver<GaussNewtonQR, corrector_t>(system,
						  std::forward<QRSolverType>(qrSolver));
}

}}
#endif  // SOLVERS_NONLINEAR_SOLVERS_CREATE_PUBLIC_API_HPP_
