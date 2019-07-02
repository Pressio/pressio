
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
class GaussNewtonQR<
  system_type, qr_solver_type,
  scalar_type, line_search_type, convergence_when_t>
  : public NonLinearSolverBase<
     GaussNewtonQR<
       system_type, qr_solver_type, scalar_type,
       line_search_type, convergence_when_t>
     >,
    public IterativeBase<scalar_type>
{
  using this_t = GaussNewtonQR<system_type, qr_solver_type,
     scalar_type, line_search_type, convergence_when_t>;

  // need to be friend of base (crpt)
  friend NonLinearSolverBase<this_t>;

  // the type of the iterative base
  using iterative_base_t = IterativeBase<scalar_type>;

  using state_t    = typename system_type::state_type;
  using residual_t = typename system_type::residual_type;
  using jacobian_t = typename system_type::jacobian_type;

  qr_solver_type qrSolver_ = {};
  residual_t res_  = {};
  jacobian_t jac_  = {};

  // to store Q^T times residual
  state_t QTResid_ = {};

  // delta is the correction
  state_t delta_     = {};

  // ytrail needed if/when line search is used
  state_t ytrial_    = {};

  // norms
  scalar_type normO_  = {};
  scalar_type normN_  = {};

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
		const state_t & yState)
    : qrSolver_{}, // default constructed
      res_(system.residual(yState)),
      jac_(system.jacobian(yState)),
      QTResid_(yState),
      delta_(yState),
      ytrial_(yState),
      normO_{0},
      normN_{0}{}

private:
  template <typename system_t>
  void solveImpl(const system_t & sys, state_t & yState){
    sys.residual(yState, res_);
    sys.jacobian(yState, jac_);

    gauss_newton_qr_solve<
      line_search_type, convergence_when_t>
      (sys, yState, ytrial_,
       res_, jac_, delta_, QTResid_,
       qrSolver_,
       iterative_base_t::maxIters_,
       iterative_base_t::tolerance_,
       normO_, normN_);
  }//end solve

};//class

}}}}//end namespace pressio::solvers::iterative::impl
#endif
