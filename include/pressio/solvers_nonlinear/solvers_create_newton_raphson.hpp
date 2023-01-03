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

#ifndef SOLVERS_NONLINEAR_SOLVERS_CREATE_NEWTON_RAPHSON_HPP_
#define SOLVERS_NONLINEAR_SOLVERS_CREATE_NEWTON_RAPHSON_HPP_

#include "solver.hpp"

namespace pressio{ namespace nonlinearsolvers{

namespace impl{
template<bool> struct _nr_choose_op;
template<> struct _nr_choose_op<true>{
  template<class ...Args> using type = impl::FusedResidualJacobianOperators<Args...>;
};
template<> struct _nr_choose_op<false>{
  template<class ...Args> using type = impl::ResidualJacobianOperators<Args...>;
};
} //end namespace impl

const std::vector<Diagnostic> defaultDiagnosticsNewtonRaphson =
  {Diagnostic::residualAbsolutel2Norm,
   Diagnostic::residualRelativel2Norm,
   Diagnostic::correctionAbsolutel2Norm,
   Diagnostic::correctionRelativel2Norm};

template<class SystemType, class LinearSolverType>
#ifdef PRESSIO_ENABLE_CXX20
requires
   (DeterminedSystemWithResidualAndJacobian<SystemType>
 || DeterminedSystemWithFusedResidualAndJacobian<SystemType>)
  && LinearSolverForNewtonRaphson<
  mpl::remove_cvref_t<LinearSolverType>,
  typename SystemType::jacobian_type,
  typename SystemType::residual_type,
  typename SystemType::state_type>
#endif
auto create_newton_raphson(const SystemType & system,
			   LinearSolverType && linSolver)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert
    (DeterminedSystemWithResidualAndJacobian<SystemType>::value
     || DeterminedSystemWithFusedResidualAndJacobian<SystemType>::value,
     "Newton-Raphson: system not satisfying the residual/jacobian concept.");

  static_assert
    (LinearSolverForNewtonRaphson<mpl::remove_cvref_t<LinearSolverType>,
     typename SystemType::jacobian_type,
     typename SystemType::residual_type,
     typename SystemType::state_type>::value,
     "Newton-Raphson: linear solver not satisfying the concept.");
#endif

  using system_type = SystemType;
  using state_t  = typename system_type::state_type;
  using r_t = typename system_type::residual_type;
  using j_t = typename system_type::jacobian_type;

  constexpr bool is_fused = DeterminedSystemWithFusedResidualAndJacobian<SystemType>;
  using op_t = typename impl::_nr_choose_op<is_fused>::template type<r_t, j_t>;
  using co_t = impl::RJCorrector<state_t, LinearSolverType>;

  op_t op(system);
  co_t co(system.createState(), std::forward<LinearSolverType>(linSolver));
  return impl::Solver<NewtonRaphson, state_t,
		      op_t, co_t>(std::move(op),
				  std::move(co),
				  defaultDiagnosticsNewtonRaphson);
}

}}
#endif  // SOLVERS_NONLINEAR_SOLVERS_CREATE_PUBLIC_API_HPP_
