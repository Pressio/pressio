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

#ifndef SOLVERS_NONLINEAR_SOLVERS_CREATE_IRLS_GAUSS_NEWTON_HPP_
#define SOLVERS_NONLINEAR_SOLVERS_CREATE_IRLS_GAUSS_NEWTON_HPP_

#include "./impl/solvers_nonlinear_compose.hpp"

namespace pressio{
namespace nonlinearsolvers{
namespace experimental{

template<class SystemType, class LinearSolverType>
#ifdef PRESSIO_ENABLE_CXX20
requires OverdeterminedSystemWithResidualAndJacobian<SystemType>
      || OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>
#endif
auto create_irls_gauss_newton(const SystemType & system,
			      LinearSolverType && linSolver)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert
    (   OverdeterminedSystemWithResidualAndJacobian<SystemType>::value
     or OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>::value,
     "Weighted Levenberg-Marquardt: system not satisfying the residual/jacobian concept.");
#endif

  using c_t = impl::ComposeIrwGaussNewton<SystemType, LinearSolverType>;
  using w_t = typename c_t::weighting_t;
  using return_t = typename c_t::type;

  w_t W(system);
  return return_t(system,
		  std::forward<LinearSolverType>(linSolver),
		  std::move(W));
}

}}}// end namespace experimental
#endif  // SOLVERS_NONLINEAR_SOLVERS_CREATE_PUBLIC_API_HPP_
