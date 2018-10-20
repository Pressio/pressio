
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_FACTORY_HPP
#define SOLVERS_EXPERIMENTAL_NONLINEAR_FACTORY_HPP

#include "../solvers_ConfigDefs.hpp"
#include "solvers_nonlinear_iterative.hpp"
#include "solvers_nonlinear_leastsquare_iterative.hpp"
#include "solvers_nonlinear_traits.hpp"
#include "solvers_policy_nonlinear_iterative.hpp"
#include "solvers_policy_nonlinear_leastsquare_iterative.hpp"
#include "solvers_policy_nonlinear_leastsquare_gauss_newton_qr.hpp"
#include "solvers_l2_vector_norm.hpp"
#include "solvers_linear_traits.hpp"
#include "../../../core/src/meta/core_meta_detection_idiom.hpp"


namespace rompp {
namespace solvers {

struct NonLinearSolvers {

  /**
   * Raise an exception while trying to create an invalid nonlinear solver.
   */
  template <
    typename NSolverT,
		typename LSolverT,
    typename std::enable_if<
      !nonlinear::details::solver_traits<NSolverT>::enabled,
      void
    >::type* = nullptr
  >
  static void createIterativeSolver() {
  	std::cerr << "Error: the nonlinear solver selected \
is not available or its name was mispelt" << std::endl;
  	assert(0);
  }
  //--------------------------------------------------------------


  /**
   * Create a nonlinear solver.
   */
  template <
    typename NSolverT,
		typename LSolverT,
    typename std::enable_if<
      nonlinear::details::solver_traits<NSolverT>::enabled,
      void
    >::type* = nullptr
  >
  static auto createIterativeSolver() {

    using solver_traits = linear::details::solver_traits<LSolverT>;

    static_assert(solver_traits::eigen_enabled && !solver_traits::direct,
		  "Error: either the linear solver is a direct one \
or is not available for linear systems defined by Eigen matrices");

    using policy_type =
      typename nonlinear::details::solver_traits<NSolverT>::solver_type;
    return NonLinearIterativeSolver<policy_type, LSolverT>();
  }
  //--------------------------------------------------------------


  /**
   * Create a nonlinear least square iterative solver
   */
  template <
    typename NSolverT,
    typename LSolverT,
    typename core::meta::enable_if_t<
      nonlinearleastsquare::details::solver_traits<NSolverT>::enabled
    >* = nullptr
  >
  static auto createNonLinearIterativeLeastSquareSolver() {

    using solver_traits = linear::details::solver_traits<LSolverT>;

    static_assert(solver_traits::eigen_enabled && !solver_traits::direct,
		  "Error: either the linear solver is a direct one \
or is not available for linear systems defined by Eigen matrices");

    using policy_type =
      typename nonlinearleastsquare::details::solver_traits<NSolverT>::solver_type;
    return NonLinearLeastSquareIterativeSolver<policy_type, LSolverT>();
  }
  //--------------------------------------------------------------

  /**
   * Create a nonlinear least square iterative solver 
   * using QR for each inner linear least square solve
   */
  template <
    typename NSolverT,
    typename core::meta::enable_if_t<
      nonlinearleastsquare::details::solver_traits<NSolverT>::enabled
    >* = nullptr
  >
  static auto createNonLinearIterativeLeastSquareQRBasedSolver(){
    using policy_type =
      typename nonlinearleastsquare::details::solver_traits<NSolverT>::solver_type;
    return NonLinearLeastSquareIterativeSolver<policy_type, void>();
  }
  //--------------------------------------------------------------

  
};

} // end namespace solvers
} // end namespace rompp
#endif
