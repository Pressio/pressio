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

#include "solver_impl.hpp"

namespace pressio{

/*To solve a determined system of nonlinear equations:
  r(x) = 0, for r \in R^n and x \in R^n.

  preconditions
  - system must bind to an object that outlives the return of this function
  - if linSolver binds an a lvalue object, it must outlive the lifetime of the object returned here

  effects
  - if linSolver binds a tempor object, it is move constructed from it
  - this function does not solve anything, only prepares for one
  - only the create* methods will be called on "sytem"
  - LinearSolver is not called, only "stored"

  post-conditions
  - the returned object own all memory allocations
*/

template<class SystemType, class LinearSolverType>
#ifdef PRESSIO_ENABLE_CXX20
  requires nonlinearsolvers::RealValuedNonlinearSystemFusingResidualAndJacobian<SystemType>
  && (Traits<typename SystemType::state_type>::rank    == 1)
  && (Traits<typename SystemType::residual_type>::rank == 1)
  && (Traits<typename SystemType::jacobian_type>::rank == 2)
  && requires(typename SystemType::state_type    & a,
	      typename SystemType::state_type    & b,
	      typename SystemType::state_type    & c,
	      typename SystemType::residual_type & r,
	      typename SystemType::jacobian_type & J,
	      nonlinearsolvers::scalar_of_t<SystemType> alpha,
	      nonlinearsolvers::scalar_of_t<SystemType> beta,
	      nonlinearsolvers::scalar_of_t<SystemType> gamma,
	      LinearSolverType && linSolver)
  {
    { ::pressio::ops::norm2(std::as_const(a)) }
      -> std::same_as< nonlinearsolvers::scalar_of_t<SystemType> >;
    { ::pressio::ops::norm2(std::as_const(r)) }
      -> std::same_as< nonlinearsolvers::scalar_of_t<SystemType> >;

    { ::pressio::ops::deep_copy(b, std::as_const(a)) };
    { ::pressio::ops::scale (a, alpha) };
    { ::pressio::ops::update(a,	alpha, std::as_const(b), beta) };
    { ::pressio::ops::update(a,	alpha, std::as_const(b), beta, std::as_const(c), gamma) };

    { linSolver.solve(std::as_const(J), std::as_const(r), a) };
  }
#endif
auto create_newton_solver(const SystemType & system,
			  LinearSolverType && linSolver)
{

  using scalar_t = nonlinearsolvers::scalar_of_t<SystemType>;
  using state_t  = typename SystemType::state_type;
  using r_t      = typename SystemType::residual_type;
  using j_t      = typename SystemType::jacobian_type;

  // A newton iteration solves: J_k delta_k = - r_k,
  // where delta_k = x_k+1 - x_k and J_k is dr_k/dx
  // so the new approximation of the solution is: x_k+1 = x_k + delta_
  using tags = std::tuple<
    nonlinearsolvers::CorrectionTag,   /*delta_k*/
    nonlinearsolvers::InitialGuessTag, /*x_0*/
    nonlinearsolvers::ResidualTag,     /*r_k*/
    nonlinearsolvers::JacobianTag,     /*J_k*/
    nonlinearsolvers::InnerSolverTag,  /*linear solver for linearization*/
    nonlinearsolvers::impl::SystemTag
    >;
  using types = std::tuple<
    state_t, state_t, r_t, j_t,
    utils::InstanceOrReferenceWrapper<LinearSolverType>,
    SystemType const *
    >;
  using registry_t = nonlinearsolvers::impl::TagBasedStaticRegistry<tags, types>;
  registry_t reg(system.createState(),
		 system.createState(),
		 system.createResidual(),
		 system.createJacobian(),
                 std::forward<LinearSolverType>(linSolver),
                 &system);

  // this is what we want to output as diagnostics during the solve
  using nonlinearsolvers::Diagnostic;
  const std::vector<Diagnostic> diagnostics =
    {Diagnostic::residualAbsolutel2Norm,
     Diagnostic::residualRelativel2Norm,
     Diagnostic::correctionAbsolutel2Norm,
     Diagnostic::correctionRelativel2Norm};

  using tag = nonlinearsolvers::impl::NewtonTag;
  return nonlinearsolvers::impl::RootFinder<tag, state_t, registry_t, scalar_t>
    (tag{}, std::move(reg), diagnostics);
}

} //end namespace pressio
#endif  // SOLVERS_NONLINEAR_SOLVERS_CREATE_PUBLIC_API_HPP_
