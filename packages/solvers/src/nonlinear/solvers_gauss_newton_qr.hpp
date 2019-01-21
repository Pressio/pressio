
#ifndef SOLVERS_GAUSS_NEWTON_QR_HPP
#define SOLVERS_GAUSS_NEWTON_QR_HPP

#include "../solvers_forward_declarations.hpp"
#include "../solvers_meta_static_checks.hpp"
#include "../base/solvers_nonlinear_base.hpp"
#include "../base/solvers_iterative_base.hpp"
#include "./impl/solvers_gauss_newton_qr_impl.hpp"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

/*
* part-specialize for when nothing about
* problem is known at compile time
*/
template <
  typename scalar_t,
  typename qr_type,
  typename line_search_t,
  typename converged_when_t
  >
class GaussNewtonQR<
  scalar_t, qr_type, line_search_t,
  converged_when_t, void, void, void, void,
  core::meta::enable_if_t<
    core::meta::is_default_constructible<qr_type>::value
    >
  > : public NonLinearSolverBase<GaussNewtonQR<scalar_t, qr_type,
					       line_search_t,
					       converged_when_t,
					       void, void,
					       void, void>>,
      public IterativeBase<scalar_t>{

  using this_t	= GaussNewtonQR<scalar_t, qr_type, line_search_t,
				converged_when_t,
				void, void, void, void>;
  using iter_base_t = IterativeBase<scalar_t>;
  using base_t  = NonLinearSolverBase<this_t>;
  friend base_t;

  scalar_t normO_ = {};
  scalar_t normN_ = {};
  qr_type qrObj	  = {};

public:
  GaussNewtonQR() = default;
  GaussNewtonQR(const GaussNewtonQR &) = delete;
  ~GaussNewtonQR() = default;

private:
  template <typename system_t>
  void solveImpl(const system_t & sys,
		 typename system_t::state_type & y){
    using state_t    = typename system_t::state_type;

    auto Resid = sys.residual(y);
    auto Jacob = sys.jacobian(y);

    state_t QTResid_(y); // hold result of Q^T residual
    QTResid_.setZero();
    state_t delta(y);
    // needed if/when doing line search: put here so that it is not
    // constructed and destroyed at each gn step
    state_t ytrial(y);

    impl::gauss_newtom_qr_solve<
      system_t, typename iter_base_t::iteration_t,
      scalar_t, qr_type, line_search_t, converged_when_t
      >(sys, y,	ytrial, Resid, Jacob, this->maxIters_,
	this->tolerance_, QTResid_, delta, qrObj,
	normO_, normN_);
  }//solve

};//class



/*
 * part-specialize for when a problem/system type is passed
 */
template <
  typename scalar_t, typename qr_type, typename line_search_t,
  typename converged_when_t, typename system_t
  >
class GaussNewtonQR<
  scalar_t, qr_type, line_search_t, converged_when_t,
  system_t, void, void, void,
  core::meta::enable_if_t<
    ::rompp::solvers::details::system_traits<system_t>::is_system and
    core::meta::is_default_constructible<qr_type>::value and
    core::meta::is_core_vector_wrapper<typename system_t::state_type>::value
    >
  > : public NonLinearSolverBase<GaussNewtonQR<scalar_t, qr_type,
					       line_search_t,
					       converged_when_t,
					       system_t, void,
					       void, void>>,
      public IterativeBase<scalar_t>
{

  using state_t    = typename system_t::state_type;
  using residual_t = typename system_t::residual_type;
  using jacobian_t = typename system_t::jacobian_type;

  using this_t	   = GaussNewtonQR<scalar_t, qr_type, line_search_t,
				   converged_when_t, system_t,
				   void, void, void>;
  using iter_base_t = IterativeBase<scalar_t>;
  using base_t	   = NonLinearSolverBase<this_t>;
  friend base_t;

  scalar_t normO_  = {};
  scalar_t normN_  = {};
  qr_type qrObj	   = {};
  state_t QTResid_ = {};
  state_t delta_   = {};
  residual_t res_  = {};
  jacobian_t jac_  = {};

  // ytrail needed if/when line search is used: put here
  // so that it is constructed only once
  state_t ytrial_  = {};

public:
  GaussNewtonQR() = delete;

  GaussNewtonQR(const system_t & system, const state_t & y)
    : QTResid_(y), delta_(y),
      res_(system.residual(y)), jac_(system.jacobian(y)),
      ytrial_(y){}

  GaussNewtonQR(const GaussNewtonQR &) = delete;
  ~GaussNewtonQR() = default;

private:
  void solveImpl(const system_t & sys, state_t & y){
    sys.residual(y, res_);
    sys.jacobian(y, jac_);

    impl::gauss_newtom_qr_solve< system_t,
				 typename iter_base_t::iteration_t,
				 scalar_t, qr_type, line_search_t,
				 converged_when_t>(sys, y, ytrial_,
						   res_, jac_,
						   this->maxIters_,
						   this->tolerance_,
						   QTResid_, delta_, qrObj,
						   normO_, normN_);
  }//end solve

};//class



/*
 * part-specialize for when the target types are passed but not system
 */
template <
  typename scalar_t, typename qr_type, typename line_search_t,
  typename converged_when_t,
  typename state_t,
  typename residual_t, typename jacobian_t
  >
class GaussNewtonQR<
  scalar_t, qr_type, line_search_t, converged_when_t,
  void, state_t, residual_t, jacobian_t,
  core::meta::enable_if_t<
    core::meta::is_default_constructible<qr_type>::value and
    core::meta::is_core_vector_wrapper<state_t>::value and
    core::meta::is_core_vector_wrapper<residual_t>::value
    >
  >
  : public NonLinearSolverBase<GaussNewtonQR<scalar_t, qr_type,
					     line_search_t,
					     converged_when_t,
					     void, state_t, residual_t,
					     jacobian_t>>,
      public IterativeBase<scalar_t>
{
  using this_t	   = GaussNewtonQR<scalar_t, qr_type, line_search_t,
				   converged_when_t, void,
				   state_t, residual_t, jacobian_t>;
  using iter_base_t = IterativeBase<scalar_t>;
  using base_t	   = NonLinearSolverBase<this_t>;
  friend base_t;

  scalar_t normO_   = {};
  scalar_t normN_   = {};
  qr_type qrObj	    = {};
  state_t QTResid_  = {};
  state_t delta_    = {};
  residual_t res_   = {};
  jacobian_t jac_   = {};

  // ytrail needed if/when line search is used: put here
  // so that it is constructed only once
  state_t ytrial_  = {};

public:
  GaussNewtonQR() = delete;

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
  GaussNewtonQR(const system_t & system, const T1 & y)
    : QTResid_(y), delta_(y),
      res_(system.residual(y)), jac_(system.jacobian(y)),
      ytrial_(y){}

  GaussNewtonQR(const GaussNewtonQR &) = delete;
  ~GaussNewtonQR() = default;

public:
  template <typename system_t>
  void solveImpl(const system_t & sys, state_t & y){
    sys.residual(y, res_);
    sys.jacobian(y, jac_);

    impl::gauss_newtom_qr_solve< system_t,
				 typename iter_base_t::iteration_t,
				 scalar_t, qr_type, line_search_t,
				 converged_when_t>(sys, y, ytrial_,
						   res_, jac_,
						   this->maxIters_,
						   this->tolerance_,
						   QTResid_, delta_,
						   qrObj, normO_, normN_);
  }//end solve

};//class


}}}}//end namespace rompp::solvers::iterative::impl
#endif
