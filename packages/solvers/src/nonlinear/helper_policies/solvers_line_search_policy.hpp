/*
//@HEADER
// ************************************************************************
//
// solvers_line_search_policy.hpp
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

#ifndef SOLVERS_IMPL_LINE_SEARCH_POLICY_HPP
#define SOLVERS_IMPL_LINE_SEARCH_POLICY_HPP

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template <typename ls_tag>
struct LineSearchHelper;

template <>
struct LineSearchHelper<gn::noLineSearch>
{
  template <typename ud_ops_t, typename scalar_t, typename ... Args>
  static void evaluate(scalar_t & alpha, Args&& ... args) {
    alpha = static_cast<scalar_t>(1);
  }
};


template <>
struct LineSearchHelper<gn::ArmijoLineSearch>{

  template <
    typename ud_ops_t,
    typename scalar_t,
    typename state_t,
    typename residual_t,
    typename jacobian_t,
    typename system_t
    >
  static void evaluate(scalar_t & alpha,
		       const state_t & y,
		       state_t & ytrial,
		       const state_t & dy,
		       residual_t & resid,
		       jacobian_t & jacob,
		       const system_t & sys)
  {
    scalar_t c1 = 1e-4;
    alpha = static_cast<scalar_t>(1);
#ifdef PRESSIO_ENABLE_DEBUG_PRINT
    ::pressio::utils::io::print_stdout("line search: Armijo rule,",
				    "c1=", c1, "\n");
#endif

    ::pressio::containers::ops::set_zero(ytrial);

    constexpr auto zero = ::pressio::utils::constants::zero<scalar_t>();
    constexpr auto one = ::pressio::utils::constants::one<scalar_t>();

    // eval obj function for current solution: f(y)
    scalar_t fy = zero;
    ComputeNormHelper::template evaluate<ud_ops_t>(resid, fy, ::pressio::solvers::Norm::L2);

    // compute J^T * Residual
    state_t jTr(y);
    ::pressio::containers::ops::set_zero(jTr);
    using jtr_prod_helper_t = JacobianTranspResProdHelper<ud_ops_t, jacobian_t>;
    // evaluate
    jtr_prod_helper_t::evaluate(jacob, resid, jTr);

    // compute dy^T J^T R (this is always a dot product)
    auto c2 = ::pressio::containers::ops::dot(dy, jTr);
    auto rhs = c1 * alpha * c2;

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
    ::pressio::utils::io::print_stdout(" f(y) =", fy, "\n");
    ::pressio::utils::io::print_stdout(" dy^T J^T R =", c2, "\n");
    ::pressio::utils::io::print_stdout(" c1*alfa*dy^T*J^T*R =", rhs, "\n");
#endif

    bool done = false;
    while (not done)
    {
#ifdef PRESSIO_ENABLE_DEBUG_PRINT
      ::pressio::utils::io::print_stdout(" backtracking: alpha =",
				      alpha, "\n");
#endif

      // update : ytrial = y + dy*alpha
      ::pressio::containers::ops::do_update(ytrial,
					    y, one,
					    dy, alpha);

      // eval function for updated step solition: f(y + alpha*dy)
      sys.residual(ytrial, resid);
      scalar_t fytrial = zero;
      ComputeNormHelper::template evaluate<ud_ops_t>(resid, fytrial, ::pressio::solvers::Norm::L2);
      auto lhs = fytrial-fy;

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
      ::pressio::utils::io::print_stdout(" f(y+alpha*dy) =", fytrial, "\n");
      ::pressio::utils::io::print_stdout(" f(y+alpha*dy)-f(y) =", lhs,
				      "; rhs =", rhs, "\n");
#endif

      // eval Armijo
      if (lhs <= rhs){
#ifdef PRESSIO_ENABLE_DEBUG_PRINT
	::pressio::utils::io::print_stdout(" lsearch done","\n");
#endif
	done = true;
      }

      // exit also when abs(fytrail-fy) < eps, leave eps = 1e-14 for now
      // change later with some machine epsilon
      if (std::abs(lhs) <= 1e-14){
#ifdef PRESSIO_ENABLE_DEBUG_PRINT
	::pressio::utils::io::print_stdout(" detected negligible",
					"change in obj f:",
					"abs(fytrail-fy) < 1e-14,",
					"exiting linsearch","\n");
#endif
	done = true;
      }

      /* convectional way to backtrack
       * this is equivalent to using beta^m instead of alpha
       * where m=0,1,2,...
       * and stopping when we find the smallest integer
       * to satisfy criterion
       */
      if (!done) alpha *= 0.5;

    }//while

#ifdef PRESSIO_ENABLE_DEBUG_PRINT
    ::pressio::utils::io::print_stdout("after line search:",
				    "alpha =", alpha, "\n");
#endif
  }//()

};


}}}} //end namespace pressio::solvers::iterative::impl
#endif
