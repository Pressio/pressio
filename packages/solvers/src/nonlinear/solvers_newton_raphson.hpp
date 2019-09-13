/*
//@HEADER
// ************************************************************************
//
// solvers_newton_raphson.hpp
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

#ifndef SOLVERS_NEWTON_RAPHSON_HPP
#define SOLVERS_NEWTON_RAPHSON_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../solvers_meta_static_checks.hpp"
#include "../../../CONTAINERS_BASIC"
#include "../../../CONTAINERS_OPS"
#include "../base/solvers_nonlinear_base.hpp"
#include "../base/solvers_iterative_base.hpp"
#include "../linear/solvers_linear_traits.hpp"

#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#endif

namespace pressio{ namespace solvers{

template <
  typename scalar_t,
  typename linear_solver_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::mpl::is_default_constructible<linear_solver_t>::value
#ifdef HAVE_PYBIND11
    or std::is_same<linear_solver_t, pybind11::object>::value
#endif
    > * = nullptr
  >
class NewtonRaphson
  : public NonLinearSolverBase<NewtonRaphson<scalar_t, linear_solver_t>>,
    public IterativeBase<NewtonRaphson<scalar_t, linear_solver_t>, scalar_t>{

  using this_t	= NewtonRaphson<scalar_t,linear_solver_t>;
  using base_t  = NonLinearSolverBase<this_t>;
  friend base_t;

  using iter_base_t  = IterativeBase<this_t, scalar_t>;
  friend iter_base_t;

  containers::default_types::uint iStep_ = {};

  // norm of the correction
  scalar_t normN_ = {};

  // residual current and previous norm
  scalar_t normRes_ = {};
  scalar_t normRes0_ = {};

  // linear solver object
  linear_solver_t linSolver_;

public:
  NewtonRaphson() = default;
  NewtonRaphson(linear_solver_t & lsO) : linSolver_{lsO}{}
  NewtonRaphson(const NewtonRaphson &) = delete;
  ~NewtonRaphson() = default;

private:
  template <typename T>
  scalar_t normOfDifference(const T & v1, const T& v2) const{
    T dVec(v1 - v2);
    return ::pressio::containers::ops::norm2(dVec);
  }

  typename iter_base_t::iteration_t getNumIterationsExecutedImpl() const {
    return iStep_;
  }


#ifdef DEBUG_PRINT
  // statically dispatch print based on whether the linear solver is iterative or not
  template <
    typename _linear_solver_t = linear_solver_t,
    ::pressio::mpl::enable_if_t<
      std::is_void<_linear_solver_t>::value == false, void
      /*::pressio::solvers::linear::details::traits<_linear_solver_t>::iterative == true*/
      > * = nullptr
    >
  void stepSummaryPrintDispatcher(const _linear_solver_t & linSolverIn,
				  const scalar_t & normRes ,
				  const scalar_t & normRes0,
				  const scalar_t & normN) const
  {
    ::pressio::utils::io::print_stdout(std::scientific,
				       "||R|| =", normRes,
				       "||R||(r) =", normRes/normRes0,
				       "||dy|| =", normN,
				       "linSolIters =", linSolverIn.getNumIterationsExecuted(),
				       "linSolFinalError =", linSolverIn.getFinalError(),
				       utils::io::reset(),
				       "\n");
  }

  template <
    typename _linear_solver_t = linear_solver_t,
    ::pressio::mpl::enable_if_t<
      ::pressio::solvers::linear::details::traits<_linear_solver_t>::direct == true
      > * = nullptr
    >
  void stepSummaryPrintDispatcher(const _linear_solver_t & linSolverIn,
  				  const scalar_t & normRes ,
  				  const scalar_t & normRes0,
  				  const scalar_t & normN) const{
    ::pressio::utils::io::print_stdout(std::scientific,
  				       "||R|| =", normRes,
  				       "||R||(r) =", normRes/normRes0,
  				       "||dy|| =", normN,
  				       utils::io::reset(),
  				       "\n");
  }
#endif

public:

  // for now, enable this for eigen only becuase it needs be cleaned
  // if we want to use it for other types
  template <
    typename system_t,
    typename state_t,
    mpl::enable_if_t<
      containers::meta::is_vector_wrapper_eigen<state_t>::value
      > * = nullptr
  >
  void solveImpl(const system_t & sys, state_t & x)
  {
#ifdef DEBUG_PRINT
    std::cout << " starting Newton-Raphson solve "
	      << " tol = " << this->tolerance_
	      << " maxIter = " << this->maxIters_
	      << std::endl;
#endif

    auto dx(x);
    state_t xOld(x);
    normN_ = {0};

    auto Residual = sys.residual(x);
    auto Jac = sys.jacobian(x);

    iStep_ = 0;

#ifdef DEBUG_PRINT
    std::cout.precision(15);
#endif
    while (++iStep_ <= this->maxIters_)
    {

#ifdef DEBUG_PRINT
      ::pressio::utils::io::print_stdout("\n");
      auto fmt = utils::io::underline();
      ::pressio::utils::io::print_stdout(fmt, "NewRaph step", iStep_,
				      utils::io::reset(), "\n");
#endif

      // solver linear system
      linSolver_.solve(Jac, Residual, dx);
      x -= dx;

      // compute norms
      normN_	= ::pressio::containers::ops::norm2(dx);
      normRes_	= ::pressio::containers::ops::norm2(Residual);
      // store initial residual norm
      if (iStep_==1) normRes0_ = normRes_;

#ifdef DEBUG_PRINT
      this->stepSummaryPrintDispatcher(linSolver_, normRes_, normRes0_, normN_);
#endif
      if (normN_ < this->tolerance_)
      	break;

      xOld = x;
      sys.residual(x, Residual);
      sys.jacobian(x, Jac);
    }
  }//solveImpl


// TODO: this needs to be cleaned
// #ifdef HAVE_PYBIND11
//   template <
//     typename system_t,
//     typename state_t,
//     mpl::enable_if_t<
//       containers::meta::is_cstyle_array_pybind11<state_t>::value and
//       mpl::is_same<system_t, pybind11::object>::value and
//       mpl::is_same<linear_solver_t, pybind11::object>::value
//       > * = nullptr
//     >
//   void solveImpl(const system_t & sys, state_t & x)
//   {
// #ifdef DEBUG_PRINT
//     std::cout << " starting Newton-Raphson solve "
// 	      << " tol = " << this->tolerance_
// 	      << " maxIter = " << this->maxIters_
// 	      << std::endl;
// #endif

//     auto dx = state_t(x.request());
//     auto xOld = state_t(dx.request());
//     auto Residual = sys.attr("residual1")(x);
//     auto Jac = sys.attr("jacobian1")(x);
//     normN_ = {0};
//     containers::default_types::uint iStep = 0;
//     std::cout.precision(15);
//     while (++iStep <= this->maxIters_)
//     {
// #ifdef DEBUG_PRINT
//       ::pressio::utils::io::print_stdout("\n");
//       auto fmt = utils::io::underline();
//       ::pressio::utils::io::print_stdout(fmt, "NewRaph step", iStep,
// 				      utils::io::reset(), "\n");
// #endif

//       xOld = x;
//       sys.attr("residual2")(x, Residual);
//       sys.attr("jacobian2")(x, Jac);

//       linSolver_.attr("solve")(Jac, dx, Residual);
//       for (auto i=0; i<x.size(); ++i)
// 	x.mutable_at(i) -= dx.at(i);

//       normN_ =::pressio::containers::ops::norm2(dx);
//       ::pressio::utils::io::print_stdout("norm(dx) =", normN_, "\n");
//       if (normN_ < this->tolerance_)
//       	break;
//     }
//   }//solveImpl
// #endif


};//class

}} //end namespace pressio::solvers
#endif
