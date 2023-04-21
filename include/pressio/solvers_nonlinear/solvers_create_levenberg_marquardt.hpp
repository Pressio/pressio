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

#ifndef SOLVERS_NONLINEAR_SOLVERS_CREATE_LEVEN_MARQ_HPP_
#define SOLVERS_NONLINEAR_SOLVERS_CREATE_LEVEN_MARQ_HPP_

#include "solvers_default_types.hpp"
#include "./impl/solvers_tagbased_registry.hpp"
#include "./impl/internal_tags.hpp"
#include "./impl/diagnostics.hpp"
#include "./impl/functions.hpp"
#include "./impl/updaters.hpp"
#include "./impl/nonlinear_least_squares.hpp"

namespace pressio{

/*
  LM minimizes the sum of squares: S(x) = (1/2) * \sum r_i(x)*r_i(x)
  where r(x) is the residual vector with r \in R^n and x \in R^k, with k < n,
  LM solves the *modified* normal equations: H delta = -g
  where:
    H = H_0 + lambda*diag(H_0), with H_0 = J^T_r*J_r
    g = J^T_r * r
    delta = x_k+1 - x_k
  The user provides a problem computing residual and jacobian
  which we can use to compute the hessian and gradient.
*/

template<class SystemType, class LinearSolverType>
#ifdef PRESSIO_ENABLE_CXX20
  requires nonlinearsolvers::RealValuedNonlinearSystemFusingResidualAndJacobian<SystemType>
  && nonlinearsolvers::valid_state_for_least_squares< typename SystemType::state_type >::value
  && (Traits<typename SystemType::state_type>::rank    == 1)
  && (Traits<typename SystemType::residual_type>::rank == 1)
  && (Traits<typename SystemType::jacobian_type>::rank == 2)
  && requires(typename SystemType::state_type & x,
	      const typename SystemType::jacobian_type & J,
	      const typename SystemType::residual_type & r,
	      nonlinearsolvers::normal_eqs_default_hessian_t<typename SystemType::state_type>  & H,
	      nonlinearsolvers::normal_eqs_default_gradient_t<typename SystemType::state_type> & g,
	      LinearSolverType && linSolver)
  {
    { ::pressio::ops::norm2(r) } -> std::same_as< nonlinearsolvers::scalar_of_t<SystemType> >;
    { ::pressio::ops::product(transpose(), nontranspose(), 1, J, 0, H) };
    { ::pressio::ops::product(transpose(), 1, J, r, 0, g) };
    { linSolver.solve(std::as_const(H), std::as_const(g), x) };
  }
#endif
auto create_levenberg_marquardt_solver(const SystemType & system,
				       LinearSolverType && linSolver)
{

  using scalar_t   = nonlinearsolvers::scalar_of_t<SystemType>;
  using state_t    = typename SystemType::state_type;
  using r_t        = typename SystemType::residual_type;
  using j_t        = typename SystemType::jacobian_type;
  using hg_default = nonlinearsolvers::normal_eqs_default_types<state_t>;
  using hessian_t  = typename hg_default::hessian_type;
  using gradient_t = typename hg_default::gradient_type;
  using lm_damp_t  = nonlinearsolvers::impl::LevenbergMarquardtDamping<scalar_t>;

  using tags = std::tuple<
    nonlinearsolvers::CorrectionTag,
    nonlinearsolvers::InitialGuessTag,
    nonlinearsolvers::ResidualTag,
    nonlinearsolvers::JacobianTag,
    nonlinearsolvers::GradientTag,
    nonlinearsolvers::HessianTag,
    nonlinearsolvers::LevenbergMarquardtUndampedHessianTag,
    nonlinearsolvers::LevenbergMarquardtDampingTag,
    nonlinearsolvers::InnerSolverTag,
    nonlinearsolvers::impl::SystemTag
    >;
  using types = std::tuple<
    state_t, state_t, r_t, j_t, gradient_t,
    hessian_t, hessian_t, lm_damp_t,
    utils::InstanceOrReferenceWrapper<LinearSolverType>,
    SystemType const *
    >;
  using registry_t = nonlinearsolvers::impl::TagBasedStaticRegistry<tags, types>;

  registry_t reg(system.createState(), system.createState(),
		 system.createResidual(), system.createJacobian(),
		 gradient_t(system.createState()),
		 hessian_t( hg_default::createHessian(system.createState()) ),
		 hessian_t( hg_default::createHessian(system.createState()) ),
		 lm_damp_t{},
                 std::forward<LinearSolverType>(linSolver),
		 &system);

  using nonlinearsolvers::Diagnostic;
  const std::vector<Diagnostic> defaultDiagnostics =
    {Diagnostic::objectiveAbsolute,
     Diagnostic::objectiveRelative,
     Diagnostic::correctionAbsolutel2Norm,
     Diagnostic::correctionRelativel2Norm,
     Diagnostic::gradientAbsolutel2Norm,
     Diagnostic::gradientRelativel2Norm};

  using tag = nonlinearsolvers::impl::LevenbergMarquardtNormalEqTag;
  return nonlinearsolvers::impl::NonLinLeastSquares<tag, state_t, registry_t, scalar_t>
    (tag{}, std::move(reg), defaultDiagnostics);
}

} // end namespace pressio
#endif  // SOLVERS_NONLINEAR_SOLVERS_CREATE_PUBLIC_API_HPP_
