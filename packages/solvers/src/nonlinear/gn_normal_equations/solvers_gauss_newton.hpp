
#ifndef SOLVERS_GAUSS_NEWTON_HPP
#define SOLVERS_GAUSS_NEWTON_HPP

#include "../../solvers_forward_declarations.hpp"
#include "../../base/solvers_nonlinear_base.hpp"
#include "../../base/solvers_iterative_base.hpp"
#include "./solvers_gauss_newton_normal_eq_impl.hpp"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template <
  typename system_type,
  typename hessian_type,
  typename linear_solver_type,
  typename scalar_type
  >
class GNHelperMixin{

protected:

  // from system_type get typedefsn
  using state_t    = typename system_type::state_type;
  using residual_t = typename system_type::residual_type;
  using jacobian_t = typename system_type::jacobian_type;

  // to compute hessian, we use a helper functor because
  // the calculation is different based on the jacobian being
  // a matrix wrapper vs a multi-vector wrapper
  using hessian_evaluator_t = HessianApproxHelper<jacobian_t>;

  GNHelperMixin() = delete;

  template <typename system_in_t,
	    typename T1 = state_t,
	    typename T2 = residual_t,
	    typename T3 = jacobian_t,
	    core::meta::enable_if_t<
	      std::is_same<T1, typename system_in_t::state_type>::value and
	      std::is_same<T2, typename system_in_t::residual_type>::value and
	      std::is_same<T3, typename system_in_t::jacobian_type>::value
	      > * = nullptr
	    >
  GNHelperMixin(const system_in_t  & system,
	       const state_t	  & yState,
	       linear_solver_type & linearSolverIn)
    : linSolver_(linearSolverIn),
      res_(system.residual(yState)),
      jac_(system.jacobian(yState)),
      hess_(hessian_evaluator_t::evaluate(jac_)),
      JTResid_(yState),
      delta_(yState),
      ytrial_(yState),
      normO_{0},
      normN_{0}{}

  linear_solver_type & linSolver_ = {};
  residual_t res_    = {};
  jacobian_t jac_    = {};
  hessian_type hess_ = {};
  state_t JTResid_   = {};

  // delta is the correction
  state_t delta_     = {};

  // ytrail needed if/when line search is used
  state_t ytrial_    = {};

  // norms
  scalar_type normO_    = {};
  scalar_type normN_    = {};
};



/* partial specialize for no observer type in templates parameters */
template <
  typename system_type,
  typename hessian_type,
  typename linear_solver_type,
  typename scalar_type,
  typename line_search_type,
  typename convergence_when_t
  >
class GaussNewton<
  system_type, hessian_type, linear_solver_type,
  scalar_type, line_search_type, convergence_when_t, void>
  : public NonLinearSolverBase<
     GaussNewton<
       system_type, hessian_type, linear_solver_type, scalar_type,
       line_search_type, convergence_when_t, void>
     >,
    public IterativeBase<scalar_type>,
    protected GNHelperMixin<system_type, hessian_type,
			   linear_solver_type, scalar_type>
{
  using this_t = GaussNewton<system_type, hessian_type, linear_solver_type,
     scalar_type, line_search_type, convergence_when_t, void>;

  // need to be friend of base (crpt)
  friend NonLinearSolverBase<this_t>;

  // the type of the iterative base
  using iterative_base_t = IterativeBase<scalar_type>;

  // mixin helper
  using gn_mixin_t = GNHelperMixin<system_type, hessian_type,
				   linear_solver_type, scalar_type>;

  using typename iterative_base_t::iteration_t;
  using typename gn_mixin_t::state_t;
  using typename gn_mixin_t::residual_t;
  using typename gn_mixin_t::jacobian_t;
  using typename gn_mixin_t::hessian_evaluator_t;

  using gn_mixin_t::linSolver_;
  using gn_mixin_t::res_;
  using gn_mixin_t::jac_;
  using gn_mixin_t::hess_;
  using gn_mixin_t::JTResid_;
  using gn_mixin_t::delta_;
  using gn_mixin_t::ytrial_;
  using gn_mixin_t::normO_;
  using gn_mixin_t::normN_;

  // dummy observer
  core::impl::empty obsObj_ = {};

public:
  GaussNewton() = delete;
  GaussNewton(const GaussNewton &) = delete;
  ~GaussNewton() = default;

  template <typename system_in_t>
  GaussNewton(const system_in_t  & system,
	      const state_t	 & yState,
	      linear_solver_type & linearSolverIn)
    :  gn_mixin_t(system, yState, linearSolverIn),
       obsObj_{}{}

private:
  template <typename system_t>
  void solveImpl(const system_t & sys, state_t & yState)
  {
    sys.residual(yState, res_);
    sys.jacobian(yState, jac_);

    gauss_newton_neq_solve<
      line_search_type, convergence_when_t>
      (sys, yState, ytrial_,
       res_, jac_, delta_, JTResid_,
       hess_, linSolver_,
       iterative_base_t::maxIters_,
       iterative_base_t::tolerance_,
       normO_, normN_, &obsObj_);

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
  typename observer_t
  >
class GaussNewton<
  system_type, hessian_type, linear_solver_type,
  scalar_type, line_search_type, convergence_when_t, observer_t>
  : public NonLinearSolverBase<
     GaussNewton<
       system_type, hessian_type, linear_solver_type, scalar_type,
       line_search_type, convergence_when_t, observer_t>
     >,
    public IterativeBase<scalar_type>,
    protected GNHelperMixin<system_type, hessian_type,
			    linear_solver_type, scalar_type>
{
  using this_t = GaussNewton<
    system_type, hessian_type, linear_solver_type, scalar_type,
    line_search_type, convergence_when_t, observer_t>;

  // need to be friend of base (crpt)
  friend NonLinearSolverBase<this_t>;

  // iterative base
  using iterative_base_t = IterativeBase<scalar_type>;

  // mixin helper
  using gn_mixin_t = GNHelperMixin<system_type, hessian_type,
				   linear_solver_type, scalar_type>;

  using typename iterative_base_t::iteration_t;
  using typename gn_mixin_t::state_t;
  using typename gn_mixin_t::residual_t;
  using typename gn_mixin_t::jacobian_t;
  using typename gn_mixin_t::hessian_evaluator_t;

  using gn_mixin_t::linSolver_;
  using gn_mixin_t::res_;
  using gn_mixin_t::jac_;
  using gn_mixin_t::hess_;
  using gn_mixin_t::JTResid_;
  using gn_mixin_t::delta_;
  using gn_mixin_t::ytrial_;
  using gn_mixin_t::normO_;
  using gn_mixin_t::normN_;

  // reference to observer object
  const observer_t & obsObj_;

public:
  GaussNewton() = delete;
  GaussNewton(const GaussNewton &) = delete;
  ~GaussNewton() = default;

  template <typename system_in_t>
  GaussNewton(const system_in_t  & system,
	      const state_t	 & yState,
	      linear_solver_type & linearSolverIn,
	      observer_t	 & obsIn)
    : gn_mixin_t(system, yState, linearSolverIn),
      obsObj_{obsIn}{}

  template <typename system_t>
  void solveImpl(const system_t & sys, state_t & yState)
  {
    sys.residual(yState, res_);
    sys.jacobian(yState, jac_);

    gauss_newton_neq_solve<
      line_search_type, convergence_when_t>
      (sys, yState, ytrial_,
       res_, jac_, delta_, JTResid_,
       hess_, linSolver_,
       iterative_base_t::maxIters_,
       iterative_base_t::tolerance_,
       normO_, normN_, &obsObj_);

  }//end solveImpl
};


}}}}//end namespace rompp::solvers::iterative::impl
#endif
