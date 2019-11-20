/*
//@HEADER
// ************************************************************************
//
// solvers_gauss_newton_qr.hpp
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

#ifndef SOLVERS_GAUSS_NEWTON_QR_HPP
#define SOLVERS_GAUSS_NEWTON_QR_HPP

#include "../../solvers_fwd.hpp"
#include "../../base/solvers_nonlinear_base.hpp"
#include "../../base/solvers_iterative_base.hpp"
#include "./solvers_gauss_newton_qr_impl.hpp"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{


/* partial specialize for no observer type in templates parameters */
template <
  typename system_type,
  typename qr_solver_type,
  typename scalar_type,
  typename line_search_type,
  typename convergence_when_t
  >
class GaussNewtonQR<system_type,
		    qr_solver_type,
		    scalar_type,
		    line_search_type,
		    convergence_when_t>
  : public NonLinearSolverBase<
     GaussNewtonQR<
       system_type, qr_solver_type, scalar_type,
       line_search_type, convergence_when_t>
     >,
    public IterativeBase<
       GaussNewtonQR<
	 system_type, qr_solver_type, scalar_type,
	 line_search_type, convergence_when_t>,
   scalar_type
   >
{
  using this_t = GaussNewtonQR<system_type, qr_solver_type,
			       scalar_type, line_search_type,
			       convergence_when_t>;

  // need to be friend of base (crpt)
  using non_lin_sol_base_t = NonLinearSolverBase<this_t>;
  friend non_lin_sol_base_t;

  // the type of the iterative base
  using iterative_base_t = IterativeBase<this_t, scalar_type>;

  using state_t    = typename system_type::state_type;
  using residual_t = typename system_type::residual_type;
  using jacobian_t = typename system_type::jacobian_type;

  qr_solver_type qrSolver_ = {};
  ::pressio::solvers::Norm normType_ = ::pressio::solvers::defaultNormType;
  residual_t residual_  = {};
  jacobian_t jacobian_  = {};
  // to store Q^T times residual
  state_t QTResid_ = {};
  // delta is the correction
  state_t correction_     = {};
  // needed if/when line search is used
  state_t trialState_    = {};

public:
  GaussNewtonQR() = delete;
  GaussNewtonQR(const GaussNewtonQR &) = delete;
  ~GaussNewtonQR() = default;

  template <typename system_in_t,
	    typename T1 = state_t,
	    typename T2 = residual_t,
	    typename T3 = jacobian_t,
	    ::pressio::mpl::enable_if_t<
	      std::is_same<T1, typename system_in_t::state_type>::value and
	      std::is_same<T2, typename system_in_t::residual_type>::value and
	      std::is_same<T3, typename system_in_t::jacobian_type>::value
	      > * = nullptr
	    >
  GaussNewtonQR(const system_in_t & system,
		const state_t & yState,
		const ::pressio::solvers::Norm normType = ::pressio::solvers::defaultNormType)
    : qrSolver_{}, // default construct
      normType_(normType),
      residual_(system.residual(yState)),
      jacobian_(system.jacobian(yState)),
      QTResid_(yState),
      correction_(yState),
      trialState_(yState){}

private:
  template <typename system_t>
  void solveImpl(const system_t & sys, state_t & yState){
    sys.residual(yState, residual_);
    sys.jacobian(yState, jacobian_);

    gauss_newton_qr_solve<
      line_search_type, convergence_when_t>
      (sys, yState, trialState_,
       residual_, jacobian_, correction_, QTResid_, qrSolver_,
       iterative_base_t::maxIters_,
       iterative_base_t::tolerance_,
       non_lin_sol_base_t::convergenceConditionDescription_, normType_);
  }//end solve

};//class

}}}}//end namespace pressio::solvers::iterative::impl
#endif
