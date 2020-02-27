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

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#ifndef SOLVERS_PY_GAUSS_NEWTON_HPP
#define SOLVERS_PY_GAUSS_NEWTON_HPP

#include "../helper_policies/solvers_converged_criterior_policy.hpp"
#include "../helper_policies/solvers_hessian_helper_policy.hpp"
#include "../helper_policies/solvers_jacob_res_product_policy.hpp"
#include "../helper_policies/solvers_norm_helper_policy.hpp"
#include "../helper_policies/solvers_line_search_policy.hpp"
#include "../helper_policies/solvers_get_matrix_size_helper.hpp"

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
  typename ops_t,
  typename when_converged_t
  >
class PyGaussNewton<
  system_t, state_t, residual_t, jacobian_t,
  hessian_t, linear_solver_t, scalar_t, ops_t, when_converged_t,
  mpl::enable_if_t<
    ::pressio::containers::meta::is_vector_wrapper_pybind<state_t>::value and
    ::pressio::containers::meta::is_vector_wrapper_pybind<residual_t>::value and
    ::pressio::containers::meta::is_matrix_wrapper_pybind<jacobian_t>::value and
    ::pressio::containers::meta::is_matrix_wrapper_pybind<hessian_t>::value
    >
  >
  : public NonLinearSolverBase<
  PyGaussNewton<
    system_t, state_t, residual_t, jacobian_t,
    hessian_t, linear_solver_t, scalar_t, ops_t, when_converged_t
    >
  >,
    public IterativeBase<
  PyGaussNewton<
    system_t, state_t, residual_t, jacobian_t,
    hessian_t, linear_solver_t, scalar_t, ops_t, when_converged_t
    >, scalar_t>
{

  static_assert(::pressio::containers::meta::is_fstyle_matrix_wrapper_pybind<jacobian_t>::value and
		::pressio::containers::meta::is_fstyle_matrix_wrapper_pybind<hessian_t>::value,
		"PyGaussNewton currently only supports jacobians and hessians \
with col-major layout. This is because we use blas and lapack for some operations.");

  static_assert(std::is_void<ops_t>::value,
		"PyGaussNewton currently only supports ops_t = void\
which means pression handles calls to compute hessian, project residual, etc. internally. \
Support for custom ops will be added later.");

  using this_t = PyGaussNewton<system_t, state_t, residual_t, jacobian_t,
			       hessian_t, linear_solver_t, scalar_t, ops_t,
			       when_converged_t>;

  // need to friend base (crpt) to grant it access to this private methods
  friend NonLinearSolverBase<this_t>;

  // the type of the iterative base
  using iterative_base_t = IterativeBase<this_t, scalar_t>;
  // the type to represent the iteration number
  using typename iterative_base_t::iteration_t;

  const ::pressio::solvers::Norm normType = ::pressio::solvers::defaultNormType;

  // policy to checking convergence
  using is_converged_t   = impl::IsConvergedHelper<when_converged_t>;

  /* policy for computing line search factor (alpha) such that
   * the update is done with y = y + alpha dy
   * alpha = 1 default when user does not want line search
   */
  using line_search_tag  = ::pressio::solvers::iterative::gn::noLineSearch;
  using lsearch_helper_t = impl::LineSearchHelper<line_search_tag>;

private:

  // here we do this conditional type because it seems when ops_t= pybind11::object
  // it only works if we copy the object. Need to figure out if we can leave ptr in all cases.
  typename std::conditional<
    mpl::is_same<ops_t, pybind11::object>::value, ops_t,
    void *
    >::type customOps_ = {};

  pybind11::object spyblas  = pybind11::module::import("scipy.linalg.blas");

  linear_solver_t linSolver_ = {};
  residual_t res_	     = {};
  jacobian_t jac_	     = {};
  hessian_t hess_	     = {};
  state_t JTR_		     = {};

  // delta is the correction
  state_t dy_		     = {};

  // ytrial needed if/when line search is used
  state_t ytrial_	     = {};

  // norms
  scalar_t normO_	     = {};
  scalar_t norm_dy_	     = {};

public:
  PyGaussNewton() = delete;
  PyGaussNewton(const PyGaussNewton &) = delete;
  PyGaussNewton(PyGaussNewton &&) = delete;
  ~PyGaussNewton() = default;

  template <
    typename _ops_t = ops_t,
    typename _system_t = system_t,
    typename _jacobian_t = jacobian_t,
    mpl::enable_if_t<
      ::pressio::containers::meta::is_matrix_wrapper_pybind<_jacobian_t>::value
      and std::is_void<_ops_t>::value
      > * = nullptr
    >
  PyGaussNewton(const _system_t	 & system,
		const typename containers::details::traits<state_t>::wrapped_t & yState,
		linear_solver_t linearSolverIn)
    : linSolver_(linearSolverIn),
      res_{system.residual(state_t(yState))},
      jac_{system.jacobian(state_t(yState))},
      hess_{jac_.extent(1), jac_.extent(1)},
      JTR_{yState.size()},
      dy_{yState.size()},
      ytrial_{yState.size()},
      normO_{0}, norm_dy_{0}
  {}

private:
  void solveImpl(const system_t & sys, state_t & y)
  {
    sys.residual(y, res_);
    sys.jacobian(y, jac_);

    constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();
    constexpr auto one = ::pressio::utils::constants::one<scalar_t>();
    constexpr auto negOne = ::pressio::utils::constants::negOne<scalar_t>();

    // alpha for taking steps
    scalar_t alpha = one;
    // storing residaul norm
    scalar_t normRes = {};
    scalar_t normRes0 = {};
    // storing projected residual norm
    scalar_t normJTRes = {};
    scalar_t normJTRes0 = {};

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
    auto ss = std::cout.precision();
    std::cout.precision(14);
    auto reset = utils::io::reset();
    auto fmt1 = utils::io::cyan() + utils::io::underline();
    const auto convString = std::string(is_converged_t::description_);
    ::pressio::utils::io::print_stdout(fmt1,
				       "PyGN normal eqns:", "criterion:",
				       convString, reset, "\n");
#endif

    // compute the initial norm of y (the state)
    impl::ComputeNormHelper::template evaluate<ops_t>(y, normO_, normType);
    norm_dy_ = {0};

    iteration_t iStep = 0;
    while (++iStep <= iterative_base_t::maxIters_)
    {

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
      ::pressio::utils::io::print_stdout("\n");
      auto fmt = utils::io::underline();
      ::pressio::utils::io::print_stdout(fmt, "PyGN step", iStep, utils::io::reset(), "\n");
#endif

      // compute norm of residual
      impl::ComputeNormHelper::template evaluate<void>(res_, normRes, normType);
      // store initial residual norm
      if (iStep==1) normRes0 = normRes;

      // //print residual and jabobian
      //pybind11::print("PyGN: residual", *res_.data(), "\n");
      // pybind11::print("PyGN: jacobian", *jac_.data(), "\n");

      //--------------------------------------------------------------
      // hessian: J^T*J
      //--------------------------------------------------------------
      constexpr auto one  = ::pressio::utils::constants::one<scalar_t>();
      constexpr auto no   = ::pressio::utils::constants::zero<int>();
      constexpr auto yes  = ::pressio::utils::constants::one<int>();
      constexpr auto transA = yes;
      constexpr auto transB = no;
      constexpr auto ovw    = yes;
      spyblas.attr("dgemm")(one, *jac_.data(), *jac_.data(),
			    zero, *hess_.data(), transA, transB, ovw);

      // // print hessian
      // pybind11::print("hessian\n", *hess_.data(), "\n");

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
      auto fmt1 = utils::io::magenta() + utils::io::bold();
      ::pressio::utils::io::print_stdout(fmt1, "GN_JSize =",
					 jac_.extent(0), jac_.extent(1), "\n");
      ::pressio::utils::io::print_stdout(fmt1, "GN_HessianSize =",
					 hess_.extent(0), hess_.extent(1),
					 utils::io::reset(), "\n");
#endif

      //--------------------------------------------------------------
      // compute gradient: J^T*res
      //--------------------------------------------------------------
      constexpr auto izero = ::pressio::utils::constants::zero<int>();
      constexpr auto ione  = ::pressio::utils::constants::one<int>();
      constexpr auto transJ= yes;
      constexpr auto ow	   = yes;
      // use -1 here to scale the result becuase of the sign convention we use
      spyblas.attr("dgemv")(negOne, *jac_.data(), *res_.data(), zero,
			    *JTR_.data(), izero, ione, izero, ione, transJ, ow);

      // // print J^T*res
      // pybind11::print("JT R");
      // pybind11::print(*JTR_.data(), "\n");

      // Norm of gradient
      impl::ComputeNormHelper::template evaluate<void>(JTR_, normJTRes, normType);
      // store the initial norm
      if (iStep==1) normJTRes0 = normJTRes;

      //--------------------------------------------------------------
      // solve normal equations
      //--------------------------------------------------------------
      linSolver_.attr("solve")(*hess_.data(), *JTR_.data(), *dy_.data());

      // compute norm of the correction
      impl::ComputeNormHelper::template evaluate<void>(dy_, norm_dy_, normType);

      // print correction
      // std::cout << "Correction dy \n" << std::endl;
      // customOps_.attr("myprint")(dy_);

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
    ::pressio::utils::io::print_stdout(std::scientific,
				    "||R|| =", normRes,
				    "||R||(r) =", normRes/normRes0,
				    "||J^T R|| =", normJTRes,
				    "||J^T R||(r) =", normJTRes/normJTRes0,
				    "||dy|| =", norm_dy_,
				    utils::io::reset(),
				    "\n");
#endif

      // //--------------------------------------------------------------
      // // update the state
      // //--------------------------------------------------------------
      // // compute multiplicative factor if needed
      // //lsearch_helper_t::evaluate(alpha, y, ytrial_, dy_, res_, jac_, sys);

      // solution update: y = y + alpha*dy;
      ::pressio::ops::do_update(y, one, dy_, alpha);

      //--------------------------------------------------------------
      // check for convergence
      //--------------------------------------------------------------
      // check convergence (whatever method user decided)
      const auto flag = is_converged_t::evaluate(y, dy_, norm_dy_, normRes, normRes0,
						 normJTRes, normJTRes0, iStep,
						 iterative_base_t::maxIters_,
						 iterative_base_t::tolerance_);

      // if we have converged, exit
      if (flag) { break; }

      // store new norm into old variable
      normO_ = norm_dy_;

      // compute residual and jacobian
      sys.residual(y, res_);
      sys.jacobian(y, jac_);
   }//loop

#if defined PRESSIO_ENABLE_DEBUG_PRINT
    std::cout.precision(ss);
    ::pressio::utils::io::print_stdout(std::fixed);
#endif

  }//end solveImpl
};

}}}//end namespace pressio::solvers::iterative
#endif
#endif






  // template <
  //   typename _ops_t = ops_t,
  //   typename _system_t = system_t,
  //   typename _jacobian_t = jacobian_t,
  //   mpl::enable_if_t<
  //     ::pressio::containers::meta::is_array_pybind11<_jacobian_t>::value and
  //     mpl::is_same<_ops_t, pybind11::object>::value
  //     > * = nullptr
  //   >
  // PyGaussNewton(const _system_t	 & system,
  // 		const state_t	 & yState,
  // 		linear_solver_t & linearSolverIn,
  // 		pybind11::object ops)
  //   : customOps_{ops},
  //     linSolver_(linearSolverIn),
  //     res_(system.residual(yState)),
  //     jac_(system.jacobian(yState)),
  //     // the hessian is square with number of rows = cols of jac
  //     // so we can initialize from just the shape
  //     // rather than calling the ops object:
  //     // hess_(customOps_.attr("multiply")(jac_, true, jac_, false)),
  //     hess_({ jac_.shape(1), jac_.shape(1) }),
  //     JTR_{ state_t(const_cast<state_t &>(yState).request()) },
  //     dy_{ state_t(const_cast<state_t &>(yState).request()) },
  //     ytrial_{ state_t(const_cast<state_t &>(yState).request()) },
  //     normO_{0}, norm_dy_{0},
  //     obsObj_{}
  // {}
