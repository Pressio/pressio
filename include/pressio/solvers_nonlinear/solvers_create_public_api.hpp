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

#ifndef SOLVERS_NONLINEAR_SOLVERS_CREATE_PUBLIC_API_HPP_
#define SOLVERS_NONLINEAR_SOLVERS_CREATE_PUBLIC_API_HPP_

#include "./impl/solvers_nonlinear_compose.hpp"

namespace pressio{ namespace nonlinearsolvers{

/*
  below we use static asserts to check constraints but this is
  not fully correct because constraints should have an impact on the
  overload resolution read this:
    https://timsong-cpp.github.io/cppwp/n4861/structure#footnote-154

  Since we cannot yet use c++20 concepts, we should enforce these
  constraints via e.g. SFINAE but that would yield bad error messages.
  So for now we decide to use static asserts to have readable error messages.
  Another point that kind of justifies this, for now, is that
  we don't really have many overloads of a given function to enable/disable
  via concepts or sfinae. So we could interpret the checks on the SystemType
  as a "mandate" rather than a constraint, so the static assert is ok.
  All this might change later if/when we decide to e.g. rename:
     create_gauss_newtonQR -> create_gauss_newton
  in such case, we would really need sfinae/concepts for the args
  to distinguish between linear solver vs qr solver.
*/

template<class SystemType, class LinearSolverType>
auto create_newton_raphson(const SystemType & system,
			   LinearSolverType && linSolver)
{
  static_assert
    (DeterminedSystemWithResidualAndJacobian<SystemType>::value or
     DeterminedSystemWithFusedResidualAndJacobian<SystemType>::value,
     "Newton-Raphson: system not satisfying the residual/jacobian concept.");

  return impl::ComposeNewtonRaphson_t<SystemType, LinearSolverType>
    (system, std::forward<LinearSolverType>(linSolver));
}

template<class SystemType, class LinearSolverType>
auto create_gauss_newton(const SystemType & system,
			 LinearSolverType && linSolver)
{

  static_assert
    (   OverdeterminedSystemWithResidualAndJacobian<SystemType>::value
     or OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>::value
     or SystemWithHessianAndGradient<SystemType>::value
     or SystemWithFusedHessianAndGradient<SystemType>::value,
     "Gauss-Newton: system does not satisfy the residual/jacobian, or hessian/gradient concepts.");

  return impl::ComposeGaussNewton_t<SystemType, LinearSolverType>
    (system, std::forward<LinearSolverType>(linSolver));
}

template<class SystemType, class LinearSolverType, class WeightingOpType>
auto create_gauss_newton(const SystemType & system,
			 LinearSolverType && linSolver,
                         WeightingOpType && weightOperator)
{

  static_assert
    (   OverdeterminedSystemWithResidualAndJacobian<SystemType>::value
     or OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>::value,
     "Weighted Gauss-Newton: system not satisfying the residual/jacobian concept.");

  return impl::ComposeGaussNewton_t<SystemType, LinearSolverType, WeightingOpType>
    (system, std::forward<LinearSolverType>(linSolver),
     std::forward<WeightingOpType>(weightOperator));
}

template<class SystemType, class QRSolverType>
auto create_gauss_newtonQR(const SystemType & system,
			   QRSolverType && qrSolver)
{

  static_assert
    (   OverdeterminedSystemWithResidualAndJacobian<SystemType>::value
     or OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>::value,
     "Gauss-Newton with QR: system not satisfying the residual/jacobian concept.");

  return impl::ComposeGaussNewtonQR_t<SystemType, QRSolverType>
    (system, std::forward<QRSolverType>(qrSolver));
}

template<class SystemType, class LinearSolverType>
auto create_levenberg_marquardt(const SystemType & system,
				LinearSolverType && linSolver)
{

  static_assert
    (   OverdeterminedSystemWithResidualAndJacobian<SystemType>::value
     or OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>::value
     or SystemWithHessianAndGradient<SystemType>::value
     or SystemWithFusedHessianAndGradient<SystemType>::value,
     "Levenberg-Marquardt: system not satisfying any of the residual/jacobian or the hessian/gradient concepts.");

  return impl::ComposeLevenbergMarquardt_t<SystemType, LinearSolverType>
    (system, std::forward<LinearSolverType>(linSolver));
}

template<class SystemType, class LinearSolverType, class WeightingOpType>
auto create_levenberg_marquardt(const SystemType & system,
				LinearSolverType && linSolver,
				WeightingOpType && weightOperator)
{

  static_assert
    (   OverdeterminedSystemWithResidualAndJacobian<SystemType>::value
     or OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>::value,
     "Weighted Levenberg-Marquardt: system not satisfying the residual/jacobian concept.");

  return impl::ComposeLevenbergMarquardt_t<SystemType, LinearSolverType, WeightingOpType>
    (system, std::forward<LinearSolverType>(linSolver),
     std::forward<WeightingOpType>(weightOperator));
}

namespace experimental{
template<class SystemType, class LinearSolverType>
auto create_irls_gauss_newton(const SystemType & system,
			      LinearSolverType && linSolver)
{

  static_assert
    (   OverdeterminedSystemWithResidualAndJacobian<SystemType>::value
     or OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>::value,
     "Weighted Levenberg-Marquardt: system not satisfying the residual/jacobian concept.");

  using c_t = impl::ComposeIrwGaussNewton<SystemType, LinearSolverType>;
  using w_t = typename c_t::weighting_t;
  using return_t = typename c_t::type;

  w_t W(system);
  return return_t(system,
		  std::forward<LinearSolverType>(linSolver),
		  std::move(W));
}
}// end namespace experimental

}}
#endif  // SOLVERS_NONLINEAR_SOLVERS_CREATE_PUBLIC_API_HPP_
