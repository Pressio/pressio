
#ifndef SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_FACTORY_HPP_
#define SOLVERS_EXPERIMENTAL_NONLINEAR_ITERATIVE_FACTORY_HPP_

#include <iostream>

#include <type_traits>
#include <Eigen/Sparse>

#include "solvers_traits.hpp"
#include "solvers_nonlinear_iterative_newton_raphson.hpp"

#include "matrix/core_matrix_traits_exp.hpp"
#include "vector/core_vector_traits_exp.hpp"


namespace solvers {

struct NonlinearIterativeSolvers {
 

  /**
   * @brief  createSolver
   * @param  An object representing a nonlinear system to be solved
   * @return A nonlinear iterative solver
   *
   * @section DESCRIPTION
   *
   * Create a nonlinear iterative solver of the specified type
   */
  template <typename SolverT>
  static auto createSolver()
  {
    typedef typename nonlinear::details::solver_traits<SolverT> solver_traits;
    typedef typename solver_traits::solver_type solver_type;
    static_assert(solver_traits::solver_exists, "Error: the nonlinear solver selected is not available");

    solver_type solver;
    return solver;
  }

};

} // end namespace solvers

#endif
