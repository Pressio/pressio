
#ifndef SOLVERS_GAUSS_NEWTON_HPP
#define SOLVERS_GAUSS_NEWTON_HPP

#include "../solvers_forward_declarations.hpp"
#include "../solvers_meta_static_checks.hpp"
#include "../base/solvers_nonlinear_base.hpp"
#include "../base/solvers_iterative_base.hpp"
#include "./impl/solvers_gauss_newton_normal_eq_impl.hpp"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

/*
* part-specialize when nothing about problem known at compile time
*/
template <
  typename scalar_t,
  typename lin_solver_tag,
  template <typename, typename> class lin_solver_t,
  typename line_search_t,
  typename converged_when_t,
  typename observer_t
  >
class GaussNewton<
  scalar_t, lin_solver_tag, lin_solver_t, line_search_t,
  converged_when_t, void, void, void, void, void, observer_t>
  : public NonLinearSolverBase<GaussNewton<scalar_t, lin_solver_tag,
					   lin_solver_t,
					   line_search_t,
					   converged_when_t,
					   void, void,
					   void, void, void,
					   observer_t>>,
      public IterativeBase<scalar_t>{

  using this_t	= GaussNewton<scalar_t, lin_solver_tag, lin_solver_t,
			      line_search_t, converged_when_t,
			      void, void, void, void, void, observer_t>;
  using iter_base_t = IterativeBase<scalar_t>;
  using base_t  = NonLinearSolverBase<this_t>;
  friend base_t;

  scalar_t normO_		= {};
  scalar_t normN_		= {};
  const observer_t * const observer_ = nullptr;

public:
  GaussNewton() = default;

  template <
    typename T = observer_t,
    core::meta::enable_if_t<
      !std::is_void<T>::value
      > * = nullptr
    >
  GaussNewton(const T & obs)
    : observer_(&obs){}

  GaussNewton(const GaussNewton &) = delete;

  ~GaussNewton() = default;

private:
  template <typename system_t>
  void solveImpl(const system_t & sys,
		 typename system_t::state_type & y){

    // note: some types are only known here because the system is passed here
    using state_t = typename system_t::state_type;
    using res_t	  = typename system_t::residual_type;
    using jac_t   = typename system_t::jacobian_type;

    static_assert(std::is_void<observer_t>::value or
		  has_method_observe_residual_when_solver_converged<observer_t, res_t>::value or
		  has_method_observe_residual_each_step<observer_t, res_t>::value,
		  "invalid observer, likely not having correct methods");

    // since I know the type of jacobian, I can declare the
    // functor that knows how to approximate the hessian: J^T*J
    using hessian_evaluator_t = HessianApproxHelper<jac_t>;

    // compure residual and jacobian
    res_t Resid = sys.residual(y);
    jac_t Jacob = sys.jacobian(y);

    // define the function and compute the hessian approx
    hessian_evaluator_t hessEvaluator;
    auto H = hessEvaluator(Jacob);

    using hessian_t = decltype(H);
    using solver_t = lin_solver_t<lin_solver_tag, hessian_t>;
    solver_t linSolver_ = {};
    static_assert(core::meta::is_default_constructible<solver_t>::value,"");

    // define variable storing J^T*res
    state_t JTR(y);

    // this is the GN correction
    state_t delta(y);

    // needed if/when doing line search: put here so that it is not
    // constructed and destroyed at each gn step
    state_t ytrial(y);

    impl::gauss_newtom_neq_solve<
      system_t, hessian_t, typename iter_base_t::iteration_t,
      scalar_t, solver_t, line_search_t,
      converged_when_t, observer_t
      >(sys, y,	ytrial, Resid, Jacob, H, JTR,
	this->maxIters_, this->tolerance_,
	delta, linSolver_, normO_, normN_, observer_);

  }//solveImpl

};//class




/*
 * part-specialize when system type is passed
 */
template <
  typename scalar_t, typename lin_solver_tag,
  template <typename,typename> class lin_solver_t,
  typename line_search_t, typename converged_when_t,
  typename system_t, typename hessian_t, typename observer_t
  >
class GaussNewton<
  scalar_t, lin_solver_tag, lin_solver_t, line_search_t,
  converged_when_t, system_t, void, void, void,
  hessian_t, observer_t,
  core::meta::enable_if_t<
    core::meta::is_core_vector_wrapper<typename system_t::state_type>::value and
    core::meta::is_core_vector_wrapper<typename system_t::residual_type>::value
    and
    (core::meta::is_core_matrix_wrapper<typename system_t::jacobian_type>::value or
     core::meta::is_core_multi_vector_wrapper<typename system_t::jacobian_type>::value)
    and
    (core::meta::is_core_matrix_wrapper<hessian_t>::value or
     core::meta::is_core_multi_vector_wrapper<hessian_t>::value)
    >
  >
  : public NonLinearSolverBase<GaussNewton<scalar_t, lin_solver_tag,
					   lin_solver_t,
					   line_search_t,
					   converged_when_t,
					   system_t, void, void,
					   void, hessian_t,
					   observer_t>>,
      public IterativeBase<scalar_t>
{

  using state_t    = typename system_t::state_type;
  using residual_t = typename system_t::residual_type;
  using jacobian_t = typename system_t::jacobian_type;
  using solverT   = lin_solver_t<lin_solver_tag, hessian_t>;

  static_assert(std::is_void<observer_t>::value or
		has_method_observe_residual_when_solver_converged<observer_t, residual_t>::value or
		has_method_observe_residual_each_step<observer_t, residual_t>::value,
		"invalid residual observer, likely not having correct methods");

  using this_t	   = GaussNewton<scalar_t, lin_solver_tag, lin_solver_t,
				 line_search_t, converged_when_t, system_t,
				 void, void, void, hessian_t, observer_t>;
  using iter_base_t = IterativeBase<scalar_t>;
  using base_t	   = NonLinearSolverBase<this_t>;
  friend base_t;

  solverT linSolver_ = {};
  scalar_t normO_    = {};
  scalar_t normN_    = {};
  state_t JTResid_   = {};
  state_t delta_     = {};
  residual_t res_    = {};
  jacobian_t jac_    = {};
  hessian_t hes_     = {};

  // ytrail needed if/when line search is used: put here
  // so that it is constructed only once
  state_t ytrial_  = {};

  const observer_t * const observer_ = nullptr;

public:
  GaussNewton() = delete;

  GaussNewton(const system_t & system, const state_t & y)
    : JTResid_(y), delta_(y),
      res_(system.residual(y)), jac_(system.jacobian(y)),
      hes_(HessianApproxHelper<jacobian_t>()(jac_)),
      ytrial_(y){}

  template <typename T = observer_t,
  	    core::meta::enable_if_t<
  	      !std::is_void<T>::value
  	      > * = nullptr
  	    >
  GaussNewton(const system_t & system,
  	      const state_t & y,
  	      const T & obs)
    : JTResid_(y), delta_(y),
      res_(system.residual(y)), jac_(system.jacobian(y)),
      hes_(HessianApproxHelper<jacobian_t>()(jac_)),
      ytrial_(y), observer_{&obs}{}

  GaussNewton(const GaussNewton &) = delete;

  ~GaussNewton() = default;

public:
  void solveImpl(const system_t & sys, state_t & y){
    sys.residual(y, res_);
    sys.jacobian(y, jac_);

    impl::gauss_newtom_neq_solve<
      system_t, hessian_t,
      typename iter_base_t::iteration_t,
      scalar_t, solverT, line_search_t,
      converged_when_t, observer_t>(sys, y, ytrial_,
    				    res_, jac_, hes_,
    				    JTResid_,
    				    this->maxIters_,
    				    this->tolerance_,
    				    delta_, linSolver_,
    				    normO_, normN_,
    				    observer_);
  }//end solve

};//class




/*
 * part-specialize for when the target types are passed but not system
 */
template <
  typename scalar_t, typename lin_solver_tag,
  template <typename,typename> class lin_solver_t,
  typename line_search_t, typename converged_when_t,
  typename state_t, typename residual_t,
  typename jacobian_t, typename hessian_t,
  typename observer_t
  >
class GaussNewton<
  scalar_t, lin_solver_tag, lin_solver_t, line_search_t,
  converged_when_t, void, state_t, residual_t,
  jacobian_t, hessian_t, observer_t,
  core::meta::enable_if_t<
    core::meta::is_core_vector_wrapper<state_t>::value and
    core::meta::is_core_vector_wrapper<residual_t>::value
    and
    (core::meta::is_core_matrix_wrapper<jacobian_t>::value or
     core::meta::is_core_multi_vector_wrapper<jacobian_t>::value)
    and
    (core::meta::is_core_matrix_wrapper<hessian_t>::value or
     core::meta::is_core_multi_vector_wrapper<hessian_t>::value)
    >
  >
  : public NonLinearSolverBase<GaussNewton<scalar_t, lin_solver_tag,
					   lin_solver_t,
					   line_search_t,
					   converged_when_t,
					   void, state_t, residual_t,
					   jacobian_t, hessian_t,
					   observer_t>>,
      public IterativeBase<scalar_t>
{

  static_assert(std::is_void<observer_t>::value or
		has_method_observe_residual_when_solver_converged<observer_t, residual_t>::value or
		has_method_observe_residual_each_step<observer_t, residual_t>::value,
		"invalid observer, likely not having correct methods");

  using solverT   = lin_solver_t<lin_solver_tag, hessian_t>;

  using this_t	   = GaussNewton<scalar_t, lin_solver_tag, lin_solver_t,
				 line_search_t, converged_when_t, void,
				 state_t, residual_t, jacobian_t,
				 hessian_t, observer_t>;
  using iter_base_t = IterativeBase<scalar_t>;
  using base_t	   = NonLinearSolverBase<this_t>;
  friend base_t;

  solverT linSolver_ = {};
  scalar_t normO_    = {};
  scalar_t normN_    = {};
  state_t JTResid_   = {};
  state_t delta_     = {};
  residual_t res_    = {};
  jacobian_t jac_    = {};
  hessian_t hes_     = {};

  // ytrail needed if/when line search is used: put here
  // so that it is constructed only once
  state_t ytrial_  = {};

  const observer_t * const observer_ = nullptr;

public:
  GaussNewton() = delete;

  template <typename system_t,
	    typename T1 = state_t,
	    typename T2 = residual_t,
	    typename T3 = jacobian_t,
	    core::meta::enable_if_t<
	      std::is_same<T1, typename system_t::state_type>::value and
	      std::is_same<T2, typename system_t::residual_type>::value and
	      std::is_same<T3, typename system_t::jacobian_type>::value
	      > * = nullptr
	    >
  GaussNewton(const system_t & system,
	      const T1 & y)
    : JTResid_(y), delta_(y),
      res_(system.residual(y)), jac_(system.jacobian(y)),
      hes_(HessianApproxHelper<jacobian_t>()(jac_)),
      ytrial_(y){}

  template <typename system_t,
  	    typename T1 = state_t,
  	    typename T2 = residual_t,
  	    typename T3 = jacobian_t,
  	    typename T4 = observer_t,
  	    core::meta::enable_if_t<
  	      std::is_same<T1, typename system_t::state_type>::value and
  	      std::is_same<T2, typename system_t::residual_type>::value and
  	      std::is_same<T3, typename system_t::jacobian_type>::value and
  	      !std::is_void<T4>::value
  	      > * = nullptr
  	    >
  GaussNewton(const system_t & system,
  	      const T1 & y,
  	      const T4 & obs)
    : JTResid_(y), delta_(y),
      res_(system.residual(y)), jac_(system.jacobian(y)),
      hes_(HessianApproxHelper<jacobian_t>()(jac_)),
      ytrial_(y), observer_{&obs}
  {}

  GaussNewton(const GaussNewton &) = delete;

  ~GaussNewton() = default;

public:
  template <typename system_t>
  void solveImpl(const system_t & sys, state_t & y){
    sys.residual(y, res_);
    sys.jacobian(y, jac_);

    impl::gauss_newtom_neq_solve
      < system_t, hessian_t,
	typename iter_base_t::iteration_t,
	scalar_t, solverT, line_search_t,
	converged_when_t, observer_t
	>(sys, y, ytrial_,
	  res_, jac_, hes_,
	  JTResid_,
	  this->maxIters_,
	  this->tolerance_,
	  delta_, linSolver_,
	  normO_, normN_, observer_);

  }//end solve

};//class


}}}}//end namespace rompp::solvers::iterative::impl
#endif
