
#ifndef SOLVERS_LINEAR_ITERATIVE_FACTORY_HPP_
#define SOLVERS_LINEAR_ITERATIVE_FACTORY_HPP_

#include <iostream>

#include <type_traits>
#include <Eigen/Sparse>

#include "solvers_traits.hpp"
#include "solvers_linear_iterative_eigen.hpp"
#include "solvers_linear_iterative_trilinos.hpp"
#include "solvers_linear_iterative_policies_trilinos.hpp"

#include "matrix/core_matrix_traits_exp.hpp"


namespace solvers {


struct LinearIterativeSolvers {
 

  template <typename SolverT,
    typename MatrixT,
    typename std::enable_if<!(
      core::details::matrix_traits<MatrixT>::matrix_class == core::details::WrappedClass::Eigen || 
      core::details::matrix_traits<MatrixT>::matrix_class == core::details::WrappedClass::Trilinos), 
      MatrixT
    >::type* = nullptr
  > 
  static void createSolver(
    MatrixT const& A
  ) {

    // Linear system solver cannot be created from given input matrix
    std::cerr << "No linear solver available for the matrix used to specify the linear system" << std::endl;
    assert(false);
  }


  template <typename SolverT,
    typename PrecT,
    typename MatrixT,
    typename std::enable_if<!(
      core::details::matrix_traits<MatrixT>::matrix_class == core::details::WrappedClass::Eigen || 
      core::details::matrix_traits<MatrixT>::matrix_class == core::details::WrappedClass::Trilinos), 
      MatrixT
    >::type* = nullptr
  >
  static void createSolver(
    MatrixT const& A
  ) {
 
    // Linear system solver cannot be created from given input matrix
    std::cerr << "No linear solver available for the matrix used to specify the linear system" << std::endl;
    assert(false);
  } 

 
  /**
   * @brief  createSolver
   * @param  A EPetra matrix object representing a linear system
   * @return An iterative linear solver for Epetra matrices based on AztecOO
   *
   * @section DESCRIPTION
   *
   * Create an instance of TrilinosLinearIterativeSolver with specified solver  
   */
  template <typename SolverT,
    typename MatrixT,
    typename std::enable_if<
      core::details::matrix_traits<MatrixT>::matrix_class == core::details::WrappedClass::Trilinos, 
      MatrixT
    >::type* = nullptr
  >
  static auto createSolver(
    MatrixT const& A
  ) {
    return LinearIterativeSolvers::template createSolver<SolverT, typename trilinos_policies::DefaultPreconditioner>(A);
  }


  /**
   * @brief  createSolver
   * @param  A EPetra matrix object representing a linear system
   * @return An iterative linear solver for Epetra matrices based on AztecOO
   *
   * @section DESCRIPTION
   *
   * Create an instance of TrilinosLinearIterativeSolver with specified solver 
   * and preconditioner types
   */
  template <typename SolverT, 
    typename PrecT,
    typename MatrixT,
    typename std::enable_if<
      core::details::matrix_traits<MatrixT>::matrix_class == core::details::WrappedClass::Trilinos, 
      MatrixT
    >::type* = nullptr
  >
  static auto createSolver(
    MatrixT const& A
  ) {

    typedef linear::details::solver_traits<SolverT> solver_traits;
    typedef linear::details::preconditioner_traits<PrecT> preconditioner_traits;

    static_assert(solver_traits::trilinos_enabled, "Solver not available for linear systems defined by Epetra matrices");
    static_assert(preconditioner_traits::trilinos_enabled, "Preconditioner not available for linear systems defined by Epetra matrices");

    TrilinosLinearIterativeSolver<typename solver_traits::trilinos_policy, 
      typename preconditioner_traits::trilinos_policy
    > solver (A);
    
    return solver; 
  }


  /**
   * @brief  createSolver
   * @param  A Eigen::SparseMatrix object representing a linear system
   * @return An Iterative linear solver for Eigen sparse matrices
   *
   * @section DESCRIPTION
   *
   * Create an instance of the appropriate sparse linear iterative solver 
   * for linear systems defined by Eigen sparse matrices
   */
  template <typename SolverT,
    typename PrecT,
    typename MatrixT,
    typename std::enable_if<
      core::details::matrix_traits<MatrixT>::matrix_class == core::details::WrappedClass::Eigen, 
      MatrixT
    >::type* = nullptr
  >
  static auto createSolver(
    MatrixT const& A
  ) {

    typedef linear::details::preconditioner_traits<PrecT> preconditioner_traits;
    static_assert(preconditioner_traits::eigen_enabled, "Preconditioner not available for linear systems defined by Eigen matrices");

    return LinearIterativeSolvers::template createSolver<SolverT>(A);
  }


  /**
   * @brief  createSolver
   * @param  A Eigen::SparseMatrix object representing a linear system
   * @return An Iterative linear solver for Eigen sparse matrices
   *
   * @section DESCRIPTION
   *
   * Create an instance of the appropriate sparse linear iterative solver 
   * for linear systems defined by Eigen sparse matrices
   */
  template <typename SolverT,
    typename MatrixT,
    typename std::enable_if<
      core::details::matrix_traits<MatrixT>::matrix_class == core::details::WrappedClass::Eigen, 
      MatrixT
    >::type* = nullptr
  >
  static auto createSolver(
    MatrixT const& A
  ) {

    typedef linear::details::solver_traits<SolverT> solver_traits;
    typedef core::details::matrix_traits<MatrixT> matrix_traits;
  
    static_assert(matrix_traits::is_sparse, "Iterative solvers in Eigen can only be used with sparse matrices");
    static_assert(solver_traits::eigen_enabled, "Solver not available for linear systems defined by Eigen matrices");
    
    typedef typename matrix_traits::wrapped_type wrapped_type;
    typedef typename solver_traits::template eigen_solver_type<wrapped_type> eigen_solver_type;
 
    EigenLinearIterativeSolver<eigen_solver_type> solver(A);
    return solver;    
  }

};

} // end namespace solvers

#endif
