
#ifndef SOLVERS_NONLINEAR_ITERATIVE_FACTORY_HPP_
#define SOLVERS_NONLINEAR_ITERATIVE_FACTORY_HPP_

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
   * @param  An object that does not represent a nonlinear system 
   * @return Fail at compile time
   */ 
  template <typename SolverT,
    typename SystemT,
    typename std::enable_if<
      !details::system_traits<SystemT>::is_system,
      SystemT
    >::type* = nullptr
  > 
  static void createSolver(
    SystemT const& A
  ) {

    // The system object used to initialize the solver is not valid
    std::cerr << "Error: the system object supplied is not valid" << std::endl;
    assert(false);
  }


  /**
   * @brief  createSolver
   * @param  An object representing a nonlinear system to be solved
   * @return A nonlinear iterative solver
   *
   * @section DESCRIPTION
   *
   * Create a nonlinear iterative solver of the specified type
   */
  template <typename SolverT,
    typename SystemT,
    typename std::enable_if<
      details::system_traits<SystemT>::is_system,
      SystemT
    >::type* = nullptr
  >
  static auto createSolver(
    SystemT const& A
  ) {

    typedef typename nonlinear::details::solver_traits<SolverT> solver_traits;
    typedef typename solver_traits::template solver_type<SystemT> solver_type;
    static_assert(solver_traits::solver_exists, "Error: the nonlinear solver selected is not available");

    solver_type solver(A);
    return solver;
  }

};

} // end namespace solvers

#endif
