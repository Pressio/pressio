
#ifndef SOLVERS_EXPERIMENTAL_LINEAR_FACTORY_HPP_
#define SOLVERS_EXPERIMENTAL_LINEAR_FACTORY_HPP_

#include <iostream>

#include <type_traits>
#include <Eigen/Sparse>

#include <solvers_traits.hpp>
#include <solvers_linear_iterative_eigen.hpp>


namespace solvers {


struct LinearSolvers {
 
  // Unsopported matrix type
  template <typename SolverT,
    typename MatrixT
  > 
  static auto createIterativeSolver(
    MatrixT const& A,
    typename std::enable_if<!(details::matrix_traits<MatrixT>::is_eigen && details::matrix_traits<MatrixT>::is_trilinos), int> = 0
  ) {

    static_assert(0, "No linear solver available for the matrix used to specify the linear system");
    return 0;

  }

  // Unsupported matrix type 
  template <typename SolverT,
    typename PrecT,
    typename MatrixT
  >
  static auto createIterativeSolver(
    MatrixT const& A,
    typename std::enable_if<!(details::matrix_traits<MatrixT>::is_eigen && details::matrix_traits<MatrixT>::is_trilinos), int> = 0
  ) {

    static_assert(0, "No linear solver available for the matrix used to specify the linear system");
    return 0;

  } 

 
  /**
   * @brief createIterativeSolver
   * @param A EPetra matrix object representing a linear system
   * @return solver Iterative linear solver for Epetra matrices based on AztecOO
   *
   * Create an instance of TrilinosLinearIterativeSolver with specified solver  
   */
  template <typename SolverT,
    typename MatrixT
  >
  static auto createIterativeSolver(
    MatrixT const& A,
    typename std::enable_if<details::matrix_traits<MatrixT>::is_trilinos, int> = 0
  ) {

    static_assert(details::solver_traits<SolverT>::trilinos_enabled, "Solver not available for linear systems defined by Eetra matrices");
    
    typedef typename core::details::matrix_traits<MatrixT>::wrapped_t wrapped_t;
  
    TrilinosLinearIterativeSolver<SolverT> solver(A);
    return solver;
  
  }


  /**
   * @brief createIterativeSolver
   * @param A EPetra matrix object representing a linear system
   * @return solver Iterative linear solver for Epetra matrices based on AztecOO
   *
   * Create an instance of TrilinosLinearIterativeSolver with specified solver 
   * and preconditioner types
   */
  template <typename SolverT, 
    typename PrecT,
    typename MatrixT
  >
  static auto createIterativeSolver(
    MatrixT const& A,
    typename std::enable_if<details::matrix_traits<MatrixT>::is_trilinos, int> = 0
  ) {

    static_assert(details::solver_traits<SolverT>::trilinos_enabled, "Solver not available for linear systems defined by Epetra matrices");
    static_assert(details::preconditioner_traits<PrecT>::trilinos_enabled, "Preconditioner not available for linear systems defined by Epetra matrices");

    typedef typename core:details::matrix_traits<MatrixT>::wrapped_t wrapped_t;

    TrilinosLinearIterativeSolver<SolvType, PrecType> solver(A);
    return solver; 

  }


  /**
   * @brief createIterativeSolver
   * @param A Eigen::SparseMatrix object representing a linear system
   * @return solver Iterative linear solver for Eigen sparse matrices
   *
   * Create an instance of the appropriate sparse linear iterative solver 
   * for linear systems defined by Eigen sparse matrices
   */
  template <typename SolverT,
    typename MatrixT
  >
  static auto createIterativeSolver( 
    MatrixT const& A,
    typename std::enable_if<details::matrix_traits<MatrixT>::is_eigen, int> = 0
  ) {

    static_assert(details::matrix_traits<MatrixT>::is_sparse, "Iterative solvers in Eigen can only be used with sparse matrices");
    static_assert(details::solver_traits<SolverT>::eigen_enabled, "Solver not available for linear systems defined by Eigen matrices");
    
    typedef typename core::details::matrix_traits<MatrixT>::wrapped_t wrapped_t;
    typedef typename details::solver_traits<SolverT>::template eigen_solver_type<wrapped_t> EigenSolverT;
 
    EigenLinearIterativeSolver<EigenSolverT> solver(A);
    return solver;    

  } 
 
};

} // end namespace solvers

#endif
