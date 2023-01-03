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

#include "solver.hpp"

namespace pressio{ namespace nonlinearsolvers{

const std::vector<Diagnostic> defaultDiagnosticsGaussNewton =
  {Diagnostic::residualAbsolutel2Norm,
   Diagnostic::residualRelativel2Norm,
   Diagnostic::correctionAbsolutel2Norm,
   Diagnostic::correctionRelativel2Norm,
   Diagnostic::gradientAbsolutel2Norm,
   Diagnostic::gradientRelativel2Norm};

namespace impl{
template<bool> struct _gn_choose_op;
template<> struct _gn_choose_op<true>{
  template<class ...Args> using type = impl::FusedResidualJacobianOperators<Args...>;
};
template<> struct _gn_choose_op<false>{
  template<class ...Args> using type = impl::ResidualJacobianOperators<Args...>;
};

template<bool> struct _gn_choose_hg_op;
template<> struct _gn_choose_hg_op<true>{
  template<class ...Args> using type = impl::FusedHessianGradientOperators<Args...>;
};
template<> struct _gn_choose_hg_op<false>{
  template<class ...Args> using type = impl::HessianGradientOperators<Args...>;
};
} //end namespace impl

template<class T, class = void>
struct NormalEquationsDefaultOperatorsTypes{
  using hessian_type  = void;
  using gradient_type = void;
};

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template<class T>
struct NormalEquationsDefaultOperatorsTypes<
  T, mpl::enable_if_t<::pressio::is_vector_eigen<T>::value> >
{
  using hessian_type = Eigen::Matrix<typename Traits<T>::scalar_type, -1, -1>;
  using gradient_type = T;
};
#endif

// ----------------------------------------------------------------
// Gauss-Newton
// ----------------------------------------------------------------
template<class SystemType, class LinearSolverType>
#ifdef PRESSIO_ENABLE_CXX20
requires
     (OverdeterminedSystemWithResidualAndJacobian<SystemType>
   || OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>)
  && LinearSolverForNonlinearLeastSquares<
       mpl::remove_cvref_t<LinearSolverType>,
       typename SystemType::state_type>
#endif
auto create_gauss_newton(const SystemType & system,
			 LinearSolverType && linSolver)
{
  using linear_solver_type = mpl::remove_cvref_t<LinearSolverType>;

#if not defined PRESSIO_ENABLE_CXX20
  static_assert
    (   OverdeterminedSystemWithResidualAndJacobian<SystemType>::value
     or OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>::value,
     "Gauss-Newton: system does not satisfy the residual/jacobian concepts.");

  static_assert(LinearSolverForNonlinearLeastSquares<
		linear_solver_type,
		typename SystemType::state_type>::value,
		"Gauss-newton: linear solver not satisfying the concept.");
#endif

  using system_type = SystemType;
  using state_t    = typename system_type::state_type;
  using r_t        = typename system_type::residual_type;
  using j_t        = typename system_type::jacobian_type;
  using hessian_t  = typename linear_solver_type::matrix_type;
  using gradient_t = state_t;

  constexpr bool is_fused = OverdeterminedSystemWithFusedResidualAndJacobian<SystemType>;
  using rj_t = typename impl::_gn_choose_op<is_fused>::template type<r_t, j_t>;
  using op_t = impl::HessianGradientDecorator<hessian_t, gradient_t, rj_t>;
  using co_t = impl::HessianGradientCorrector<state_t, LinearSolverType>;

  op_t op(system);
  co_t co(system.createState(), std::forward<LinearSolverType>(linSolver));
  return impl::Solver<
    GaussNewton, state_t, op_t, co_t>(std::move(op),
				      std::move(co),
				      defaultDiagnosticsGaussNewton);
}

template<class SystemType, class LinearSolverType>
#ifdef PRESSIO_ENABLE_CXX20
requires
  (SystemWithHessianAndGradient<SystemType>
  || SystemWithFusedHessianAndGradient<SystemType>)
  && LinearSolverForNonlinearLeastSquares<
       mpl::remove_cvref_t<LinearSolverType>,
       typename SystemType::state_type>
#endif
auto create_gauss_newton(const SystemType & system,
			 LinearSolverType && linSolver)
{
  using linear_solver_type = mpl::remove_cvref_t<LinearSolverType>;

#if not defined PRESSIO_ENABLE_CXX20
  static_assert
    (   SystemWithHessianAndGradient<SystemType>::value
     or SystemWithFusedHessianAndGradient<SystemType>::value,
     "Gauss-Newton: system does not satisfy the hessian/gradient concepts.");
  static_assert(LinearSolverForNonlinearLeastSquares<
		linear_solver_type,
		typename SystemType::state_type>::value,
		"Gauss-newton: linear solver not satisfying the concept.");
#endif

  using system_type = SystemType;
  using state_t     = typename system_type::state_type;
  using hessian_t   = typename system_type::hessian_type;
  using gradient_t  = typename system_type::gradient_type;
  using norm_type   = typename system_type::residual_norm_type;

  constexpr bool is_fused = SystemWithFusedHessianAndGradient<SystemType>;
  using op_t = typename impl::_gn_choose_hg_op<is_fused>::template type<hessian_t, gradient_t, norm_type>;
  using co_t = impl::HessianGradientCorrector<state_t, LinearSolverType>;

  op_t op(system);
  co_t co(system.createState(), std::forward<LinearSolverType>(linSolver));
  return impl::Solver<
    GaussNewton, state_t, op_t, co_t>(std::move(op),
				      std::move(co),
				      defaultDiagnosticsGaussNewton);
}

}}
#endif  // SOLVERS_NONLINEAR_SOLVERS_CREATE_PUBLIC_API_HPP_
