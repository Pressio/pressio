/*
//@HEADER
// ************************************************************************
//
// solvers_lm_neq_res_jac_api.hpp
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

#ifndef SOLVERS_LM_NORMAL_EQ_RES_JAC_API_HPP
#define SOLVERS_LM_NORMAL_EQ_RES_JAC_API_HPP

#include "./solvers_lm_neq_solve_impl.hpp"
#include "../../../helpers/solvers_norm_dispatcher.hpp"
#include "../../../helpers/solvers_hessian_dispatcher.hpp"
#include "../../../helpers/solvers_gradient_dispatcher.hpp"

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template <
  typename system_type,
  typename hessian_type,
  typename linear_solver_type,
  typename scalar_type,
  typename ud_ops_t
  >
class LMHelperMixin
{

protected:
  /* --- type aliasing--- */
  using state_t    = typename system_type::state_type;
  using residual_t = typename system_type::residual_type;
  using jacobian_t = typename system_type::jacobian_type;

  /* --- members --- */
  linear_solver_type & linSolver_ = {};
  ::pressio::solvers::Norm normType_ = ::pressio::solvers::defaultNormType;
  ::pressio::solvers::iterative::impl::HessianDispatcher<ud_ops_t> hessianDispatcher_ = {};
  ::pressio::solvers::iterative::impl::GradientDispatcher<ud_ops_t> gradientDispatcher_ = {};
  ::pressio::solvers::iterative::impl::NormDispatcher<ud_ops_t> normDispatcher_ = {};
  residual_t residual_    = {};
  jacobian_t jacobian_    = {};
  hessian_type hessian_	  = {};
  state_t gradient_	  = {};
  state_t correction_     = {};
  state_t trialState_    = {};

protected:
  LMHelperMixin() = delete;
  ~LMHelperMixin() = default;

  template <
    typename system_in_t,
    typename T1 = state_t,
    typename T2 = residual_t,
    typename T3 = jacobian_t,
    typename = ::pressio::mpl::enable_if_t<
      std::is_same<T1, typename system_in_t::state_type>::value and
      std::is_same<T2, typename system_in_t::residual_type>::value and
      std::is_same<T3, typename system_in_t::jacobian_type>::value and
      std::is_void<ud_ops_t>::value
      >
    >
  LMHelperMixin(const system_in_t  & system,
		const state_t	  & yState,
		linear_solver_type & linearSolverIn,
		const ::pressio::solvers::Norm normType)
    : linSolver_(linearSolverIn),
      normType_(normType),
      residual_(system.residual(yState)),
      jacobian_(system.jacobian(yState)),
      hessian_(hessianDispatcher_.template evaluate<jacobian_t, hessian_type>(jacobian_)),
      gradient_(yState),
      correction_(yState),
      trialState_(yState)
  {}

  template <
    typename system_in_t,
    typename T1 = state_t,
    typename T2 = residual_t,
    typename T3 = jacobian_t,
    typename _ud_ops_t = ud_ops_t,
    typename = ::pressio::mpl::enable_if_t<
      std::is_same<T1, typename system_in_t::state_type>::value and
      std::is_same<T2, typename system_in_t::residual_type>::value and
      std::is_same<T3, typename system_in_t::jacobian_type>::value and
      !std::is_void<ud_ops_t>::value
      >
    >
  LMHelperMixin(const system_in_t  & system,
		const state_t	  & yState,
		linear_solver_type & linearSolverIn,
		const ::pressio::solvers::Norm normType,
		const _ud_ops_t & udOps)
    : linSolver_(linearSolverIn),
      normType_(normType),
      hessianDispatcher_{&udOps},
      gradientDispatcher_{&udOps},
      normDispatcher_{&udOps},
      residual_(system.residual(yState)),
      jacobian_(system.jacobian(yState)),
      hessian_(hessianDispatcher_.template evaluate<jacobian_t, hessian_type>(jacobian_)),
      gradient_(yState),
      correction_(yState),
      trialState_(yState)
  {}
};


template <
  typename system_t,
  typename hessian_t,
  typename linear_solver_t,
  typename scalar_t,
  typename lm_schedule_policy_tag,
  typename when_converged_t,
  typename resid_obs_t,
  typename ud_ops_t,
  typename enable = void
  >
class LMNormalEqResJacApi;


/* partial specialize for no observer type in templates parameters */
template <
  typename system_type,
  typename hessian_type,
  typename linear_solver_type,
  typename scalar_type,
  typename lm_schedule_policy_tag,
  typename convergence_when_t,
  typename ud_ops_t
  >
class LMNormalEqResJacApi<
  system_type,
  hessian_type,
  linear_solver_type,
  scalar_type,
  lm_schedule_policy_tag,
  convergence_when_t,
  void,
  ud_ops_t
  >
  : public NonLinearSolverBase<
     LMNormalEqResJacApi<system_type, hessian_type, linear_solver_type, scalar_type,
				  lm_schedule_policy_tag, convergence_when_t, void, ud_ops_t>
     >,
    public IterativeBase< LMNormalEqResJacApi<system_type, hessian_type,
						       linear_solver_type, scalar_type,
						       lm_schedule_policy_tag, convergence_when_t,
						       void, ud_ops_t>, scalar_type>,
    protected LMHelperMixin<system_type, hessian_type, linear_solver_type, scalar_type, ud_ops_t>
{

  using this_t = LMNormalEqResJacApi<system_type, hessian_type, linear_solver_type,
					      scalar_type, lm_schedule_policy_tag, convergence_when_t,
					      void, ud_ops_t>;

  // need to be friend of base (crpt)
  using non_lin_sol_base_t = NonLinearSolverBase<this_t>;
  friend non_lin_sol_base_t;

  // the iterative base type
  using iterative_base_t = IterativeBase<this_t, scalar_type>;

  // mixin helper
  using lm_mixin_t = LMHelperMixin<system_type, hessian_type, linear_solver_type, scalar_type, ud_ops_t>;

  using typename iterative_base_t::iteration_t;
  using typename lm_mixin_t::state_t;
  using typename lm_mixin_t::residual_t;
  using typename lm_mixin_t::jacobian_t;

  using lm_mixin_t::linSolver_;
  using lm_mixin_t::normType_;
  using lm_mixin_t::residual_;
  using lm_mixin_t::jacobian_;
  using lm_mixin_t::hessian_;
  using lm_mixin_t::gradient_;
  using lm_mixin_t::correction_;
  using lm_mixin_t::trialState_;
  using lm_mixin_t::hessianDispatcher_;
  using lm_mixin_t::gradientDispatcher_;
  using lm_mixin_t::normDispatcher_;

  // dummy observer
  utils::impl::empty obsObj_ = {};

  using lm_schedule_policy_t = pressio::solvers::iterative::impl::LMSchedule<lm_schedule_policy_tag,scalar_type>;
  lm_schedule_policy_t LMSchedule_;

public:
  LMNormalEqResJacApi() = delete;
  LMNormalEqResJacApi(const LMNormalEqResJacApi &) = delete;
  ~LMNormalEqResJacApi() = default;

  template <
    typename system_in_t, typename _ud_ops_t = ud_ops_t,
    typename = mpl::enable_if_t< std::is_void<_ud_ops_t>::value >
    >
  LMNormalEqResJacApi(const system_in_t  & system,
			       const state_t	 & yState,
			       linear_solver_type & linearSolverIn,
			       const ::pressio::solvers::Norm normType = ::pressio::solvers::defaultNormType)
    :  lm_mixin_t(system, yState, linearSolverIn, normType),
       obsObj_{},
       LMSchedule_(){}

  template <
    typename system_in_t, typename _ud_ops_t = ud_ops_t,
    typename = mpl::enable_if_t< std::is_void<_ud_ops_t>::value >,
    typename lm_schedule_t
    >
  LMNormalEqResJacApi(const system_in_t  & system,
			       const state_t	 & yState,
			       linear_solver_type & linearSolverIn,
             lm_schedule_t & lmSchedule,
			       const ::pressio::solvers::Norm normType = ::pressio::solvers::defaultNormType)
    :  lm_mixin_t(system, yState, linearSolverIn, normType),
       obsObj_{},
       LMSchedule_(lmSchedule){}


  template <
    typename system_in_t, typename _ud_ops_t = ud_ops_t,
    typename = mpl::enable_if_t< !std::is_void<_ud_ops_t>::value >
    >
  LMNormalEqResJacApi(const system_in_t  & system,
			       const state_t	 & yState,
			       linear_solver_type & linearSolverIn,
			       const _ud_ops_t & udOps,
			       const ::pressio::solvers::Norm normType = ::pressio::solvers::defaultNormType)
    :  lm_mixin_t(system, yState, linearSolverIn, normType, udOps),
       obsObj_{},
       LMSchedule_(){}

  template <
    typename system_in_t, typename _ud_ops_t = ud_ops_t,
    typename = mpl::enable_if_t< !std::is_void<_ud_ops_t>::value >,
    typename lm_schedule_t
    >
  LMNormalEqResJacApi(const system_in_t  & system,
			       const state_t	 & yState,
			       linear_solver_type & linearSolverIn,
			       const _ud_ops_t & udOps,
             lm_schedule_t & lmSchedule,
			       const ::pressio::solvers::Norm normType = ::pressio::solvers::defaultNormType)
    :  lm_mixin_t(system, yState, linearSolverIn, normType, udOps),
       obsObj_{},
       LMSchedule_(){}


private:
  template <typename system_t>
  void solveImpl(const system_t & sys, state_t & yState)
  {
    LMSchedule_.reset();
    sys.residual(yState, residual_);
    sys.jacobian(yState, jacobian_);

    lm_neq_solve<lm_schedule_policy_tag, convergence_when_t>
      (sys, yState, trialState_,
       residual_, jacobian_, correction_, gradient_,
       hessian_, linSolver_,
       iterative_base_t::maxIters_,
       iterative_base_t::tolerance_,
       &obsObj_,
       non_lin_sol_base_t::convergenceConditionDescription_,
       normType_,
       hessianDispatcher_, gradientDispatcher_, normDispatcher_,LMSchedule_);
  }//end solveImpl
};



/* non-void observer type is passed to templates parameters */
template <
  typename system_type,
  typename hessian_type,
  typename linear_solver_type,
  typename scalar_type,
  typename lm_schedule_policy_tag,
  typename convergence_when_t,
  typename observer_t,
  typename ud_ops_t
  >
class LMNormalEqResJacApi<
  system_type,
  hessian_type,
  linear_solver_type,
  scalar_type,
  lm_schedule_policy_tag,
  convergence_when_t,
  observer_t,
  ud_ops_t
  >
  : public NonLinearSolverBase< LMNormalEqResJacApi<system_type, hessian_type,
							     linear_solver_type, scalar_type,
							     lm_schedule_policy_tag, convergence_when_t,
							     observer_t, ud_ops_t>>,
    public IterativeBase< LMNormalEqResJacApi<system_type, hessian_type,
						       linear_solver_type, scalar_type,
						       lm_schedule_policy_tag, convergence_when_t,
						       observer_t, ud_ops_t>, scalar_type>,
    protected LMHelperMixin<system_type, hessian_type, linear_solver_type, scalar_type, ud_ops_t>
{

  using this_t = LMNormalEqResJacApi<system_type, hessian_type,
					      linear_solver_type, scalar_type,
					      lm_schedule_policy_tag, convergence_when_t,
					      observer_t, ud_ops_t>;

  // need to be friend of base (crpt)
  using non_lin_sol_base_t = NonLinearSolverBase<this_t>;
  friend non_lin_sol_base_t;

  // iterative base
  using iterative_base_t = IterativeBase<this_t, scalar_type>;

  // mixin helper
  using lm_mixin_t = LMHelperMixin<system_type, hessian_type, linear_solver_type, scalar_type, ud_ops_t>;

  using typename iterative_base_t::iteration_t;
  using typename lm_mixin_t::state_t;
  using typename lm_mixin_t::residual_t;
  using typename lm_mixin_t::jacobian_t;

  using lm_mixin_t::linSolver_;
  using lm_mixin_t::normType_;
  using lm_mixin_t::residual_;
  using lm_mixin_t::jacobian_;
  using lm_mixin_t::hessian_;
  using lm_mixin_t::gradient_;
  using lm_mixin_t::correction_;
  using lm_mixin_t::trialState_;
  using lm_mixin_t::hessianDispatcher_;
  using lm_mixin_t::gradientDispatcher_;
  using lm_mixin_t::normDispatcher_;

  // reference to observer object
  const observer_t & obsObj_;

  using lm_schedule_policy_t = pressio::solvers::iterative::impl::LMSchedule<lm_schedule_policy_tag,scalar_type>;
  lm_schedule_policy_t LMSchedule_;

public:
  LMNormalEqResJacApi() = delete;
  LMNormalEqResJacApi(const LMNormalEqResJacApi &) = delete;
  ~LMNormalEqResJacApi() = default;

  template <typename system_in_t>
  LMNormalEqResJacApi(const system_in_t  & system,
			       const state_t	 & yState,
			       linear_solver_type & linearSolverIn,
			       observer_t	 & obsIn,
			       const ::pressio::solvers::Norm normType = ::pressio::solvers::defaultNormType)
    : lm_mixin_t(system, yState, linearSolverIn, normType),
      obsObj_{obsIn},
      LMSchedule_(){}

  template <typename system_in_t, typename lm_schedule_t>
  LMNormalEqResJacApi(const system_in_t  & system,
			       const state_t	 & yState,
			       linear_solver_type & linearSolverIn,
			       observer_t	 & obsIn,
             lm_schedule_t & lmSchedule,
			       const ::pressio::solvers::Norm normType = ::pressio::solvers::defaultNormType)
    : lm_mixin_t(system, yState, linearSolverIn, normType),
      obsObj_{obsIn},
      LMSchedule_(lmSchedule){}


  template <typename system_t>
  void solveImpl(const system_t & sys, state_t & yState)
  {
    LMSchedule_.reset();
    sys.residual(yState, residual_);
    sys.jacobian(yState, jacobian_);
    lm_neq_solve<
      lm_schedule_policy_tag, convergence_when_t
      >(sys, yState, trialState_,
	residual_, jacobian_, correction_, gradient_,
	hessian_, linSolver_,
	iterative_base_t::maxIters_,
	iterative_base_t::tolerance_,
	&obsObj_,
	non_lin_sol_base_t::convergenceConditionDescription_,
	normType_,
	hessianDispatcher_, gradientDispatcher_, normDispatcher_,LMSchedule_);
  }//end solveImpl

};

}}}}//end namespace pressio::solvers::nonlinear::impl
#endif
