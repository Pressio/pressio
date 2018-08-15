
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_FACTORY_HPP
#define SOLVERS_EXPERIMENTAL_NONLINEAR_FACTORY_HPP

#include <iostream>
#include <type_traits>

#include "solvers_nonlinear_iterative.hpp"
#include "solvers_nonlinear_traits.hpp"


namespace solvers {
	
struct NonLinearSolvers {

  /**
   * @brief Raise an assertion as the nonlinear solver specified is invalid
   *
   */
  template <
    typename SolverT,
    typename std::enable_if<
      !nonlinear::details::solver_traits<SolverT>::enabled,
      int
    >::type* = nullptr
  >
  static void createSolver() {
  	std::cerr << "Error: the nonlinear solver selected is not available or its name was mispelt" << std::endl;
  	assert(0);
  }
  
  /**
   * @brief Create a valid non linear solver
   *
   */
  template <
    typename SolverT,
    typename std::enable_if<
      nonlinear::details::solver_traits<SolverT>::enabled,
      int
    >::type* = nullptr
  >
  static auto createSolver() {
  	return NonLinearIterativeSolver<void>();
  }

};

} // end namespace solvers

#endif