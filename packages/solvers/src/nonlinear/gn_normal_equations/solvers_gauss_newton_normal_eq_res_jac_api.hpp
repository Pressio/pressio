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

#ifndef SOLVERS_GAUSS_NEWTON_NORMAL_EQ_RES_JAC_API_HPP
#define SOLVERS_GAUSS_NEWTON_NORMAL_EQ_RES_JAC_API_HPP

#include "./solvers_gauss_newton_normal_eq_impl.hpp"

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template <
  typename system_type,
  typename hessian_type,
  typename linear_solver_type,
  typename scalar_type,
  typename ud_ops_t = void
  >
class GNHelperMixin
{

protected:
  /* --- type aliasing--- */
  using state_t    = typename system_type::state_type;
  using residual_t = typename system_type::residual_type;
  using jacobian_t = typename system_type::jacobian_type;
  // to compute hessian, we use a helper functor because
  // the calculation is different based on the jacobian being
  // a matrix wrapper vs a multi-vector wrapper
  using hessian_evaluator_t = HessianApproxHelper<ud_ops_t, jacobian_t>;

  /* --- members --- */
  linear_solver_type & linSolver_ = {};
  ::pressio::solvers::Norm normType_ = ::pressio::solvers::defaultNormType;
  residual_t residual_    = {};
  jacobian_t jacobian_    = {};
  hessian_type hessian_	  = {};
  state_t gradient_	  = {};
  state_t correction_     = {};
  state_t trialState_    = {};

protected:
  GNHelperMixin() = delete;
  ~GNHelperMixin() = default;

  template <
    typename system_in_t,
    typename T1 = state_t,
    typename T2 = residual_t,
    typename T3 = jacobian_t,
    ::pressio::mpl::enable_if_t<
      std::is_same<T1, typename system_in_t::state_type>::value and
      std::is_same<T2, typename system_in_t::residual_type>::value and
      std::is_same<T3, typename system_in_t::jacobian_type>::value
      > * = nullptr
    >
  GNHelperMixin(const system_in_t  & system,
		const state_t	  & yState,
		linear_solver_type & linearSolverIn,
		const ::pressio::solvers::Norm normType)
    : linSolver_(linearSolverIn),
      normType_(normType),
      residual_(system.residual(yState)),
      jacobian_(system.jacobian(yState)),
      hessian_(hessian_evaluator_t::template evaluate<hessian_type>(jacobian_)),
      gradient_(yState),
      correction_(yState),
      trialState_(yState)
  {}
};



/* partial specialize for no observer type in templates parameters */
template <
  typename system_type,
  typename hessian_type,
  typename linear_solver_type,
  typename scalar_type,
  typename line_search_type,
  typename convergence_when_t,
  typename ud_ops_t
  >
class GaussNewtonNormalEqResJacApi<
  system_type,
  hessian_type,
  linear_solver_type,
  scalar_type,
  line_search_type,
  convergence_when_t,
  void,
  ud_ops_t
  >
  : public NonLinearSolverBase<
     GaussNewtonNormalEqResJacApi<system_type, hessian_type, linear_solver_type, scalar_type,
				  line_search_type, convergence_when_t, void, ud_ops_t>
     >,
    public IterativeBase< GaussNewtonNormalEqResJacApi<system_type, hessian_type,
						       linear_solver_type, scalar_type,
						       line_search_type, convergence_when_t,
						       void, ud_ops_t>, scalar_type>,
    protected GNHelperMixin<system_type, hessian_type, linear_solver_type, scalar_type, ud_ops_t>
{

  using this_t = GaussNewtonNormalEqResJacApi<system_type, hessian_type, linear_solver_type,
					      scalar_type, line_search_type, convergence_when_t,
					      void, ud_ops_t>;

  // need to be friend of base (crpt)
  using non_lin_sol_base_t = NonLinearSolverBase<this_t>;
  friend non_lin_sol_base_t;

  // the iterative base type
  using iterative_base_t = IterativeBase<this_t, scalar_type>;

  // mixin helper
  using gn_mixin_t = GNHelperMixin<system_type, hessian_type, linear_solver_type, scalar_type, ud_ops_t>;

  using typename iterative_base_t::iteration_t;
  using typename gn_mixin_t::state_t;
  using typename gn_mixin_t::residual_t;
  using typename gn_mixin_t::jacobian_t;
  using typename gn_mixin_t::hessian_evaluator_t;

  using gn_mixin_t::linSolver_;
  using gn_mixin_t::normType_;
  using gn_mixin_t::residual_;
  using gn_mixin_t::jacobian_;
  using gn_mixin_t::hessian_;
  using gn_mixin_t::gradient_;
  using gn_mixin_t::correction_;
  using gn_mixin_t::trialState_;

  // dummy observer
  utils::impl::empty obsObj_ = {};

public:
  GaussNewtonNormalEqResJacApi() = delete;
  GaussNewtonNormalEqResJacApi(const GaussNewtonNormalEqResJacApi &) = delete;
  ~GaussNewtonNormalEqResJacApi() = default;

  template <typename system_in_t>
  GaussNewtonNormalEqResJacApi(const system_in_t  & system,
			       const state_t	 & yState,
			       linear_solver_type & linearSolverIn,
			       const ::pressio::solvers::Norm normType = ::pressio::solvers::defaultNormType)
    :  gn_mixin_t(system, yState, linearSolverIn, normType),
       obsObj_{}{}

private:
  template <typename system_t>
  void solveImpl(const system_t & sys, state_t & yState)
  {
    sys.residual(yState, residual_);
    sys.jacobian(yState, jacobian_);

    gauss_newton_neq_solve<
      ud_ops_t, line_search_type, convergence_when_t>
      (sys, yState, trialState_,
       residual_, jacobian_, correction_, gradient_,
       hessian_, linSolver_,
       iterative_base_t::maxIters_,
       iterative_base_t::tolerance_,
       &obsObj_,
       non_lin_sol_base_t::convergenceConditionDescription_,
       normType_);
  }//end solveImpl
};



/* non-void observer type is passed to templates parameters */
template <
  typename system_type,
  typename hessian_type,
  typename linear_solver_type,
  typename scalar_type,
  typename line_search_type,
  typename convergence_when_t,
  typename observer_t,
  typename ud_ops_t
  >
class GaussNewtonNormalEqResJacApi<
  system_type,
  hessian_type,
  linear_solver_type,
  scalar_type,
  line_search_type,
  convergence_when_t,
  observer_t,
  ud_ops_t
  >
  : public NonLinearSolverBase< GaussNewtonNormalEqResJacApi<system_type, hessian_type,
							     linear_solver_type, scalar_type,
							     line_search_type, convergence_when_t,
							     observer_t, ud_ops_t>>,
    public IterativeBase< GaussNewtonNormalEqResJacApi<system_type, hessian_type,
						       linear_solver_type, scalar_type,
						       line_search_type, convergence_when_t,
						       observer_t, ud_ops_t>, scalar_type>,
    protected GNHelperMixin<system_type, hessian_type, linear_solver_type, scalar_type, ud_ops_t>
{

  using this_t = GaussNewtonNormalEqResJacApi<system_type, hessian_type,
					      linear_solver_type, scalar_type,
					      line_search_type, convergence_when_t,
					      observer_t, ud_ops_t>;

  // need to be friend of base (crpt)
  using non_lin_sol_base_t = NonLinearSolverBase<this_t>;
  friend non_lin_sol_base_t;

  // iterative base
  using iterative_base_t = IterativeBase<this_t, scalar_type>;

  // mixin helper
  using gn_mixin_t = GNHelperMixin<system_type, hessian_type, linear_solver_type, scalar_type, ud_ops_t>;

  using typename iterative_base_t::iteration_t;
  using typename gn_mixin_t::state_t;
  using typename gn_mixin_t::residual_t;
  using typename gn_mixin_t::jacobian_t;
  using typename gn_mixin_t::hessian_evaluator_t;

  using gn_mixin_t::linSolver_;
  using gn_mixin_t::normType_;
  using gn_mixin_t::residual_;
  using gn_mixin_t::jacobian_;
  using gn_mixin_t::hessian_;
  using gn_mixin_t::gradient_;
  using gn_mixin_t::correction_;
  using gn_mixin_t::trialState_;

  // reference to observer object
  const observer_t & obsObj_;

public:
  GaussNewtonNormalEqResJacApi() = delete;
  GaussNewtonNormalEqResJacApi(const GaussNewtonNormalEqResJacApi &) = delete;
  ~GaussNewtonNormalEqResJacApi() = default;

  template <typename system_in_t>
  GaussNewtonNormalEqResJacApi(const system_in_t  & system,
			       const state_t	 & yState,
			       linear_solver_type & linearSolverIn,
			       observer_t	 & obsIn,
			       const ::pressio::solvers::Norm normType = ::pressio::solvers::defaultNormType)
    : gn_mixin_t(system, yState, linearSolverIn, normType),
      obsObj_{obsIn}{}

  template <typename system_t>
  void solveImpl(const system_t & sys, state_t & yState)
  {
    sys.residual(yState, residual_);
    sys.jacobian(yState, jacobian_);

    gauss_newton_neq_solve<
      ud_ops_t, line_search_type, convergence_when_t
      >(sys, yState, trialState_,
	residual_, jacobian_, correction_, gradient_,
	hessian_, linSolver_,
	iterative_base_t::maxIters_,
	iterative_base_t::tolerance_,
	&obsObj_,
	non_lin_sol_base_t::convergenceConditionDescription_,
	normType_);
  }//end solveImpl

};

}}}}//end namespace pressio::solvers::iterative::impl
#endif
