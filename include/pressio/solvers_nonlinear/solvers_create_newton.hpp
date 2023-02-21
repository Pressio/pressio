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

namespace pressio{ namespace nonlinearsolvers{

/*
  Create a solver for solving determined systems of
  nonlinear equations of the form: r(x) = 0 , for r \in R^n and x \in R^n.

  preconditions
  - system must bind to an object that outlives the result of this function
  - if linSolver binds to a lvalue object, that must also outlive
    the lifetime of the object returned here

  effects
  - if linSolver binds a tempor object, it is move constructed from it
  - not solve is actaully performed, only preparing for it
  - only the create* methods will be called on "sytem"
  - LinearSolver is not called at all

  post-conditions
  - the returned object own all memory allocations
*/
template<class SystemType, class LinearSolverType>
#ifdef PRESSIO_ENABLE_CXX20
requires
    (DeterminedRealValuedSystemWithResidualAndJacobian<SystemType>
  || DeterminedRealValuedSystemWithFusedResidualAndJacobian<SystemType>)
  && SolverForNewtonStepOf< mpl::remove_cvref_t<LinearSolverType>, SystemType >
  && requires(typename SystemType::state_type & x,
	      typename SystemType::state_type & b,
	      typename SystemType::state_type & c,
	      typename SystemType::residual_type & r,
	      scalar_of_t<SystemType> alpha)
  {
    { ::pressio::ops::norm2(x) } -> std::same_as< scalar_of_t<SystemType> >;
    { ::pressio::ops::norm2(r) } -> std::same_as< scalar_of_t<SystemType> >;
    { ::pressio::ops::deep_copy(b, x) };
    { ::pressio::ops::scale(x, scalar_of_t<SystemType>{}) };
    { ::pressio::ops::update(x, scalar_of_t<SystemType>{},
			     b, scalar_of_t<SystemType>{}) };
    { ::pressio::ops::update(x, scalar_of_t<SystemType>{},
			     b, scalar_of_t<SystemType>{},
			     c, scalar_of_t<SystemType>{}) };
  }
#endif
auto create_newton(const SystemType & system,
		   LinearSolverType && linSolver)
{

  // since the system meets the required concepts,
  // find the nested types of the various operators
  using state_t = typename SystemType::state_type;
  using r_t     = typename SystemType::residual_type;
  using j_t     = typename SystemType::jacobian_type;
  using scalar_t  = scalar_trait_t<state_t>;

  // A newton iteration solves:  J_k delta_k = - r_k
  // where: delta_k = x_k+1 - x_k and J_k is dr_k/dx
  // so the new approximation of the solution is: x_k+1 = x_k + delta_
 // the registry contains all the data we need for the solve
  using registry_t = impl::TagBasedStaticRegistry<
    std::tuple<CorrectionTag,		/*delta_k*/
	       InitialGuessTag,		/*x_0*/
	       ResidualTag,		/*r_k*/
	       JacobianTag,		/*J_k*/
	       InnerSolverTag		/*linear solver for linearization*/
	       >,
    std::tuple<state_t, state_t, r_t, j_t,
	       utils::InstanceOrReferenceWrapper<LinearSolverType>>
    >;
  registry_t reg(system.createState(),
		 system.createState(),
		 system.createResidual(),
		 system.createJacobian(),
		 std::forward<LinearSolverType>(linSolver));

  // this is what we want to output as diagnostics during the solve
  const std::vector<Diagnostic> diagnostics =
    {Diagnostic::residualAbsolutel2Norm,
     Diagnostic::residualRelativel2Norm,
     Diagnostic::correctionAbsolutel2Norm,
     Diagnostic::correctionRelativel2Norm};

  using tag = impl::NewtonTag;
  return impl::RootFinder<tag, state_t, registry_t, scalar_t>
    (tag{}, std::move(reg), diagnostics);
}

}}
#endif  // SOLVERS_NONLINEAR_SOLVERS_CREATE_PUBLIC_API_HPP_
