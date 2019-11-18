/*
//@HEADER
// ************************************************************************
//
// solvers_gauss_newton.hpp
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

#ifndef SOLVERS_GAUSS_NEWTON_HPP
#define SOLVERS_GAUSS_NEWTON_HPP

#include "../../solvers_fwd.hpp"
#include "../../base/solvers_nonlinear_base.hpp"
#include "../../base/solvers_iterative_base.hpp"
#include "./solvers_gauss_newton_normal_eq_impl.hpp"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{
namespace experimental{

template <
  typename system_type,
  typename linear_solver_type,
  typename scalar_type,
  typename line_search_type,
  typename convergence_when_t
  >
class GaussNewtonJtjJtrApi<
  system_type,
  linear_solver_type,
  scalar_type,
  line_search_type,
  convergence_when_t
  >
  : public NonLinearSolverBase</*...*/>,
    public IterativeBase</*...*/>
{
  using this_t = /*...*/;

  // need to be friend of base (crpt)
  using non_lin_sol_base_t = NonLinearSolverBase<this_t>;
  friend non_lin_sol_base_t;

  // iterative base
  using iterative_base_t = IterativeBase<this_t, scalar_type>;

  // typedefs from the system
  using state_t    = typename system_type::state_type;
  using hessian_t  = typename system_type::hessian_type;
  using proj_res_t = typename system_type::projected_residual_type;

  // --- data members ---
  linear_solver_type & linSolver_ = {};
  hessian_type hess_		  = {};
  proj_res_t JTResid_		  = {};
  // delta is the correction
  state_t delta_		  = {};
  // ytrail needed if/when line search is used
  state_t ytrial_		  = {};
  // norms
  scalar_type normO_		  = {};
  scalar_type normN_		  = {};

public:
  GaussNewtonJtjJtrApi() = delete;
  GaussNewtonJtjJtrApi(const GaussNewtonJtjJtrApi &) = delete;
  ~GaussNewtonJtjJtrApi() = default;

  GaussNewtonJtjJtrApi(const system_type  & system,
		       const state_t	  & stateIn,
		       linear_solver_type & linearSolverIn)
    : linSolver_(linearSolverIn),
      // J^T J and J^T R constructed from the system
      hess_( system.createHessianObject(stateIn) ),
      JTResid_( system.createProjectedResidualObject(stateIn) ),
      delta_(yState),
      ytrial_(yState),
      normO_{0},
      normN_{0}
  {}

private:
  void solveImpl(const system_type & sys, state_t & yState)
  {
    /* ... */
  }
};

}}}}}//end namespace pressio::solvers::iterative::impl::experimental
#endif
