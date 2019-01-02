
#ifndef SOLVERS_LEASTSQUARES_GAUSS_NEWTON_QR_HPP
#define SOLVERS_LEASTSQUARES_GAUSS_NEWTON_QR_HPP

#include "../solvers_forward_declarations.hpp"
#include "../solvers_meta_static_checks.hpp"
#include "solvers_gauss_newton_qr_impl.hpp"

namespace rompp{ namespace solvers{

/*
* part-specialize for when nothing about
* problem is known at compile time
*/
template <typename scalar_t, typename qr_type>
class GaussNewtonQR<
  scalar_t, qr_type, void, void, void, void,
  core::meta::enable_if_t<
    core::meta::is_default_constructible<qr_type>::value
    >
  >{

  scalar_t nonLinearTolerance_       = 1e-6;
  scalar_t normO_ 		     = {};
  scalar_t normN_ 		     = {};
  qr_type qrObj			     = {};
  using uint_t			     = core::default_types::uint;
  uint_t maxNonLinearIterations_     = 500;

public:
  GaussNewtonQR() = default;
  GaussNewtonQR(const GaussNewtonQR &) = delete;
  ~GaussNewtonQR() = default;

public:
  void setMaxNonLinearIterations(uint_t maxNonLinearIterations) {
    maxNonLinearIterations_ = maxNonLinearIterations;
  }

  void setNonLinearTolerance(scalar_t nonLinearTolerance) {
    nonLinearTolerance_ = std::abs(nonLinearTolerance);
  }

  template <typename system_t,
	    core::meta::enable_if_t<
	      core::meta::is_core_vector_wrapper<typename system_t::state_type>::value
	      > * =nullptr
	    >
  void solve(const system_t & sys, typename system_t::state_type & x){
    using state_t    = typename system_t::state_type;

    auto Resid = sys.residual(x);
    auto Jacob = sys.jacobian(x);

    state_t QTResid_(x); // hold result of Q^T residual
    QTResid_.setZero();
    state_t dx(x);

    ::rompp::solvers::impl::gauss_newtom_qr_solve(sys, x,
    						  Resid, Jacob,
    						  maxNonLinearIterations_,
    						  nonLinearTolerance_,
    						  QTResid_, dx, qrObj,
						  normO_, normN_);
  }//solve

};//class



/*
 * part-specialize for when a problem/system type is passed
 */
template <typename scalar_t, typename qr_type, typename system_t>
class GaussNewtonQR<
  scalar_t, qr_type, system_t, void, void, void,
  core::meta::enable_if_t<
    ::rompp::solvers::details::system_traits<system_t>::is_system and
    core::meta::is_default_constructible<qr_type>::value and
    core::meta::is_core_vector_wrapper<typename system_t::state_type>::value
    >
  >{

  using state_t    = typename system_t::state_type;
  using residual_t = typename system_t::residual_type;
  using jacobian_t = typename system_t::jacobian_type;
  using uint_t     = core::default_types::uint;

  scalar_t nonLinearTolerance_       = 1e-6;
  scalar_t normO_ 		     = {};
  scalar_t normN_ 		     = {};
  uint_t maxNonLinearIterations_     = 500;

  qr_type qrObj			     = {};
  state_t QTResid_		     = {};
  state_t delta_		     = {};
  residual_t res_		     = {};
  jacobian_t jac_		     = {};

public:
  GaussNewtonQR() = delete;

  GaussNewtonQR(const system_t & system, const state_t & y)
    : QTResid_(y), delta_(y),
      res_(system.residual(y)), jac_(system.jacobian(y)) {}

  GaussNewtonQR(const GaussNewtonQR &) = delete;
  ~GaussNewtonQR() = default;

public:
  void setMaxNonLinearIterations(uint_t maxNonLinearIterations) {
    maxNonLinearIterations_ = maxNonLinearIterations;
  }

  void setNonLinearTolerance(scalar_t nonLinearTolerance) {
    nonLinearTolerance_ = std::abs(nonLinearTolerance);
  }

  void solve(system_t & sys, state_t & x){
    sys.residual(x, res_);
    sys.jacobian(x, jac_);
    ::rompp::solvers::impl::gauss_newtom_qr_solve(sys, x,
    						  res_, jac_,
    						  maxNonLinearIterations_,
    						  nonLinearTolerance_,
    						  QTResid_, delta_, qrObj,
    						  normO_, normN_);
  }//end solve

};//class






/*
 * part-specialize for when the target types are passed but not system
 */
template <typename scalar_t, typename qr_type,
	  typename state_t, typename residual_t, typename jacobian_t>
class GaussNewtonQR<
  scalar_t, qr_type, void, state_t, residual_t, jacobian_t,
  core::meta::enable_if_t<
    core::meta::is_default_constructible<qr_type>::value and
    core::meta::is_core_vector_wrapper<state_t>::value and
    core::meta::is_core_vector_wrapper<residual_t>::value
    >
  >{

  using uint_t     = core::default_types::uint;

  scalar_t nonLinearTolerance_       = 1e-6;
  scalar_t normO_ 		     = {};
  scalar_t normN_ 		     = {};
  uint_t maxNonLinearIterations_     = 500;

  qr_type qrObj			     = {};
  state_t QTResid_		     = {};
  state_t delta_		     = {};
  residual_t res_		     = {};
  jacobian_t jac_		     = {};

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
      res_(system.residual(y)), jac_(system.jacobian(y)) {}

  GaussNewtonQR(const GaussNewtonQR &) = delete;
  ~GaussNewtonQR() = default;

public:
  void setMaxNonLinearIterations(uint_t maxNonLinearIterations) {
    maxNonLinearIterations_ = maxNonLinearIterations;
  }

  void setNonLinearTolerance(scalar_t nonLinearTolerance) {
    nonLinearTolerance_ = std::abs(nonLinearTolerance);
  }

  template <typename system_t>
  void solve(system_t & sys, state_t & x){
    sys.residual(x, res_);
    sys.jacobian(x, jac_);
    ::rompp::solvers::impl::gauss_newtom_qr_solve(sys, x,
    						  res_, jac_,
    						  maxNonLinearIterations_,
    						  nonLinearTolerance_,
    						  QTResid_, delta_, qrObj,
    						  normO_, normN_);
  }//end solve

};//class



}} //end namespace rompp::solvers
#endif
