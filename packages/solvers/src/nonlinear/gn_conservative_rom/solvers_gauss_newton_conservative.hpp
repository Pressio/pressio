
#ifndef SOLVERS_GAUSS_NEWTON_CONSERVATIVE_HPP
#define SOLVERS_GAUSS_NEWTON_CONSERVATIVE_HPP

#include "../../solvers_forward_declarations.hpp"
//#include "../../solvers_meta_static_checks.hpp"
#include "../../base/solvers_nonlinear_base.hpp"
#include "../../base/solvers_iterative_base.hpp"
#include "./solvers_gauss_newton_normal_eq_conservative_impl.hpp"

namespace rompp{ namespace solvers{ namespace iterative{ namespace hacked{

/*
 * part-specialize when system type is passed
 */
template <
  typename scalar_t, typename lin_solver_tag,
  template <typename,typename> class lin_solver_t,
  typename line_search_t, typename converged_when_t,
  typename system_t, typename cbar_t
  >
class GaussNewtonConservative<
  scalar_t, lin_solver_tag, lin_solver_t, line_search_t,
  converged_when_t, system_t, cbar_t,
  ::rompp::mpl::enable_if_t<
    core::meta::is_vector_wrapper_eigen<typename system_t::state_type>::value and
    core::meta::is_core_vector_wrapper<typename system_t::residual_type>::value
    and
    (core::meta::is_core_matrix_wrapper<typename system_t::jacobian_type>::value or
     core::meta::is_core_multi_vector_wrapper<typename system_t::jacobian_type>::value)
    and
     core::meta::is_core_multi_vector_wrapper<cbar_t>::value
    >
  >
  : public NonLinearSolverBase<
  GaussNewtonConservative<
    scalar_t, lin_solver_tag, lin_solver_t, line_search_t,
    converged_when_t, system_t, cbar_t>>,
    public IterativeBase<scalar_t>
{

  using eig_dyn_mat= Eigen::Matrix<scalar_t, -1, -1>;
  using mat_t= rompp::core::Matrix<eig_dyn_mat>;
  using state_t    = typename system_t::state_type;
  using residual_t = typename system_t::residual_type;
  using jacobian_t = typename system_t::jacobian_type;
  using solverT    = lin_solver_t<lin_solver_tag, mat_t>;

  using this_t	   = GaussNewtonConservative<scalar_t, lin_solver_tag, lin_solver_t,
					     line_search_t, converged_when_t, system_t, cbar_t>;
  using iter_base_t = IterativeBase<scalar_t>;
  using base_t	   = NonLinearSolverBase<this_t>;
  friend base_t;

  residual_t res_    = {};
  jacobian_t jac_    = {};
  residual_t cbarTlambda_ = {};
  const cbar_t & cbarT_	= {};

  solverT linSolver_ = {};
  scalar_t normO_    = {};
  scalar_t normN_    = {};
  state_t delta_     = {};
  mat_t jTj_     = {};
  mat_t jTcbarT_   = {};
  mat_t cbarJ_     = {};
  mat_t zero_	     = {};
  state_t jTr2_ = {};
  state_t cbarR_ = {};
  mat_t A_ = {};
  state_t b_ = {};
  state_t lambda_ = {};
  state_t y2_ = {};

  // ytrail needed if/when line search is used: put here
  // so that it is constructed only once
  state_t ytrial_  = {};

public:
  GaussNewtonConservative() = delete;

  /* the system passed in does not have to be the same as the
   * class template parameter BUT the two systems have to compatible
   * in the sense that their return types of residual and jacobian
   * have to be the same. maybe we should get rid of this and
   * create an overload where residual and jacobian types are
   * passed explicitly.
   */
  template <typename system_in_t = system_t>
  GaussNewtonConservative(const system_in_t & system,
			  const state_t & y,
			  const cbar_t & cbarT)
    : res_(system.residual(y)),
      jac_(system.jacobian(y)),
      cbarTlambda_(system.residual(y)),//this is just to initialize it
      cbarT_(cbarT)
  {
    // number of cols of Cbar
    const auto nlambda = cbarT_.globalNumVectors();

    jTj_.resize( jac_.globalNumVectors(), jac_.globalNumVectors() );
    jTcbarT_.resize( jac_.globalNumVectors(), nlambda );
    cbarJ_.resize( nlambda, jac_.globalNumVectors() );
    zero_.resize( nlambda, nlambda);
    zero_.setZero();

    jTr2_.resize(y.size());
    cbarR_.resize(nlambda);

    A_.resize( jac_.globalNumVectors()+nlambda, jac_.globalNumVectors()+nlambda);
    b_.resize(y.size() + nlambda);
    lambda_.resize(nlambda);
    lambda_.putScalar( static_cast<scalar_t>(0) );
    y2_.resize(y.size() + nlambda);
    ytrial_.resize(y2_.size());
    delta_.resize(y2_.size());

    y2_.data()->block(0, 0, y.size(), 1) = *y.data();
    y2_.data()->block(y.size(), 0, lambda_.size(), 1) = *lambda_.data();
  }

  GaussNewtonConservative(const GaussNewtonConservative &) = delete;
  ~GaussNewtonConservative() = default;

public:
  template <typename system_in_t = system_t>
  void solveImpl(const system_in_t & sys, state_t & y){
    sys.residual(y, res_);
    sys.jacobian(y, jac_);

    impl::gauss_newtom_neq_conserv_solve<
      system_in_t,
      typename iter_base_t::iteration_t,
      scalar_t, solverT, line_search_t,
      converged_when_t, cbar_t, mat_t
      >(sys,
	y,
	ytrial_,
	res_, jac_,
	this->maxIters_,
	this->tolerance_,
	delta_, linSolver_,
	normO_, normN_,
	cbarT_,
	jTj_, jTcbarT_,
	cbarJ_,	zero_,
	cbarTlambda_, jTr2_,
	cbarR_,	A_,
	b_, lambda_,
	y2_);

  }//end solve

};//class


}}}}//end namespace rompp::solvers::iterative::hacked
#endif
