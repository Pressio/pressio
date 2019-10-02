/*
//@HEADER
// ************************************************************************
//
// solvers_py_gauss_newton.hpp
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

#ifdef HAVE_PYBIND11
#ifndef SOLVERS_PY_GAUSS_NEWTON_HPP
#define SOLVERS_PY_GAUSS_NEWTON_HPP

#include "../../solvers_fwd.hpp"
#include "../../base/solvers_nonlinear_base.hpp"
#include "../../base/solvers_iterative_base.hpp"
#include "./solvers_gauss_newton_normal_eq_impl.hpp"

namespace pressio{ namespace solvers{ namespace iterative{

/*
 * for interfacing with python
 * for time being, no-line search
*/
template <
  typename system_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t,
  typename hessian_t,
  typename linear_solver_t,
  typename scalar_t,
  typename when_converged_t
  >
class PyGaussNewton<
  system_t, state_t, residual_t, jacobian_t,
  hessian_t, linear_solver_t, scalar_t, when_converged_t,
  mpl::enable_if_t<
    ::pressio::containers::meta::is_array_pybind11<state_t>::value and
    ::pressio::containers::meta::is_array_pybind11<residual_t>::value and
    ::pressio::containers::meta::is_array_pybind11<jacobian_t>::value and
    ::pressio::containers::meta::is_array_pybind11<hessian_t>::value
    >
  >
  : public NonLinearSolverBase<
  PyGaussNewton<
    system_t, state_t, residual_t, jacobian_t,
    hessian_t, linear_solver_t, scalar_t, when_converged_t
    >
  >,
    public IterativeBase<
  PyGaussNewton<
    system_t, state_t, residual_t, jacobian_t,
    hessian_t, linear_solver_t, scalar_t, when_converged_t
    >, scalar_t>
{
  static_assert( mpl::not_same<system_t, pybind11::object>::value,
		 "The PyGaussNewton is supposed to have state/res/jac = pybind11::array_t, but the system_t should NOT be a pybind11::object to contain a Python object. If you need to use a GaussNewton solver on numpy strucrtures, just use one ready off the shelf. This PyGaussNewton is tailored so that you can use an Application class written in Python and use the time-integrators within Pressio as well as the ROMs methods inside Pressio.");

  using this_t = PyGaussNewton<
    system_t, state_t, residual_t, jacobian_t,
    hessian_t, linear_solver_t, scalar_t, when_converged_t>;

  // need to friend base (crpt) to grant it access to this private methods
  friend NonLinearSolverBase<this_t>;

  // the type of the iterative base
  using iterative_base_t = IterativeBase<this_t, scalar_t>;

  using typename iterative_base_t::iteration_t;

  // this references a Python class where we search for all needed ops
  pybind11::object pythonOps_ = {};

  linear_solver_t linSolver_ = {};
  residual_t res_    = {};
  jacobian_t jac_    = {};
  hessian_t hess_    = {};
  state_t JTR_       = {};

  // delta is the correction
  state_t dy_     = {};

  // ytrail needed if/when line search is used
  state_t ytrial_    = {};

  // norms
  scalar_t normO_    = {};
  scalar_t norm_dy_    = {};

  // dummy observer
  utils::impl::empty obsObj_ = {};

public:
  PyGaussNewton() = delete;
  PyGaussNewton(const PyGaussNewton &) = delete;
  ~PyGaussNewton() = default;

  PyGaussNewton(const system_t	 & system,
		const state_t	 & yState,
		linear_solver_t & linearSolverIn,
		pybind11::object ops)
    : pythonOps_{ops},
      linSolver_(linearSolverIn),
      res_(system.residual(yState)),
      jac_(system.jacobian(yState)),
      hess_(pythonOps_.attr("multiply1")(jac_, jac_, true)),
      JTR_{ state_t(const_cast<state_t &>(yState).request()) },
      dy_{ state_t(const_cast<state_t &>(yState).request()) },
      ytrial_{ state_t(const_cast<state_t &>(yState).request()) },
      normO_{0}, norm_dy_{0},
      obsObj_{}
  {}

private:

  void solveImpl(const system_t & sys, state_t & y)
  {
    sys.residual(y, res_);
    sys.jacobian(y, jac_);

    // find out which norm to use
    using norm_t = typename impl::NormSelectorHelper<when_converged_t>::norm_t;

    // policy for evaluating the norm of a vector
    using norm_evaluator_t = impl::ComputeNormHelper<norm_t>;

    // policy to checking convergence
    using is_converged_t = impl::IsConvergedHelper<when_converged_t>;

    /* policy for computing line search factor (alpha) such that
     * the update is done with y = y + alpha dy
     * alpha = 1 default when user does not want line search
     */
    using line_search_tag  = ::pressio::solvers::iterative::gn::noLineSearch;
    using lsearch_helper_t = impl::LineSearchHelper<line_search_tag>;

    // alpha for taking steps
    scalar_t alpha = {};
    // storing residaul norm
    scalar_t normRes = {};
    scalar_t normRes0 = {};
    // storing projected residual norm
    scalar_t normJTRes = {};
    scalar_t normJTRes0 = {};

    constexpr auto one = ::pressio::utils::constants::one<scalar_t>();
    constexpr auto negOne = ::pressio::utils::constants::negOne<scalar_t>();

#ifdef DEBUG_PRINT
    auto ss = std::cout.precision();
    std::cout.precision(14);
    auto reset = utils::io::reset();
    auto fmt1 = utils::io::cyan() + utils::io::underline();
    const auto convString = std::string(is_converged_t::description_);
    ::pressio::utils::io::print_stdout(fmt1, "PyGN normal eqns:", "criterion:",
				    convString, reset, "\n");
#endif

    // compute the initial norm of y (the state)
    norm_evaluator_t::evaluate(y, normO_);
    norm_dy_ = {0};

    iteration_t iStep = 0;
    while (++iStep <= iterative_base_t::maxIters_)
    {
#ifdef DEBUG_PRINT
      ::pressio::utils::io::print_stdout("\n");
      auto fmt = utils::io::underline();
      ::pressio::utils::io::print_stdout(fmt, "PyGN step", iStep,
				      utils::io::reset(), "\n");
#endif

      // residual norm for current state
      norm_evaluator_t::evaluate(res_, normRes);
      // std::cout << "residual" << std::endl;
      // pythonOps_.attr("myprint")(res_);

      // store initial residual norm
      if (iStep==1) normRes0 = normRes;

      // compute LHS: J^T*J
      pythonOps_.attr("multiply2")(jac_, jac_, hess_, true);
      // std::cout << "\n";
      // std::cout << "jacobian" << std::endl;
      // pythonOps_.attr("myprint")(jac_);
      // std::cout << "\n";
#ifdef DEBUG_PRINT
      std::cout << "hessian" << std::endl;
      pythonOps_.attr("myprint")(hess_);
#endif

#ifdef DEBUG_PRINT
      auto fmt1 = utils::io::magenta() + utils::io::bold();
      ::pressio::utils::io::print_stdout(fmt1, "GN_JSize =",
				      jac_.shape()[0], jac_.shape()[1],
				      "\n");
      ::pressio::utils::io::print_stdout(fmt1, "GN_HessianSize =",
				    hess_.shape()[0], hess_.shape()[1],
				    utils::io::reset(), "\n");
#endif

      // compute RHS: J^T*res
      pythonOps_.attr("multiply2")(jac_, res_, JTR_, true);
      pythonOps_.attr("scale")(JTR_, negOne);

      // norm of projected residual
      norm_evaluator_t::evaluate(JTR_, normJTRes);
      // store initial residual norm
      if (iStep==1) normJTRes0 = normJTRes;

      // solve normal equations
      linSolver_.attr("solve")(hess_, JTR_, dy_);

      // compute norm of the correction
      norm_evaluator_t::evaluate(dy_, norm_dy_);

#ifdef DEBUG_PRINT
      ::pressio::utils::io::print_stdout(std::scientific,
				      "||R|| =", normRes,
				      "||R||(r) =", normRes/normRes0,
				      "||dy|| =", norm_dy_,
				      utils::io::reset(),
				      "\n");
#endif

      // compute multiplicative factor if needed
      lsearch_helper_t::evaluate(alpha, y, ytrial_, dy_, res_, jac_, sys);

      // solution update: y = y + alpha*dy;
      ::pressio::containers::ops::do_update(y, one, dy_, alpha);
      // std::cout << "PyGN PrintUpdateSol" << std::endl;
      // pythonOps_.attr("myprint")(dy_);
      // pythonOps_.attr("myprint")(y);
      // std::cout << std::endl;

      // check convergence (whatever method user decided)
      const auto flag = is_converged_t::evaluate(y, dy_, norm_dy_,
						 normRes, normRes0,
						 normJTRes, normJTRes0,
						 iStep,
						 iterative_base_t::maxIters_,
						 iterative_base_t::tolerance_);

      // if we have converged, exit
      if (flag) {
	break;
      }

      // store new norm into old variable
      normO_ = norm_dy_;

      // compute residual and jacobian
      sys.residual(y, res_);
      sys.jacobian(y, jac_);

      }//loop

#if defined DEBUG_PRINT
    std::cout.precision(ss);
    ::pressio::utils::io::print_stdout(std::fixed);
#endif

  }//end solveImpl
};


}}}//end namespace pressio::solvers::iterative
#endif
#endif
