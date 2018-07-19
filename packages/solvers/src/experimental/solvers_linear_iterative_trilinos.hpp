
#ifndef SOLVERS_TRILINOS_LINEAR_ITERATIVE_HPP_
#define SOLVERS_TRILINOS_LINEAR_ITERATIVE_HPP_

#include <iostream>
#include <memory>
#include <type_traits>

#include "AztecOO.h"
#include "Epetra_RowMatrix.h"

#include "matrix/core_matrix_traits_exp.hpp"
#include "vector/core_vector_traits_exp.hpp"

#include "solvers_linear_iterative_base.hpp"


namespace solvers {

// Fwd declaration
struct LinearIterativeSolvers;


/**
 * @brief Class that implements LinearIterativeSolverBase for Trilinos
 */
template<typename SolverT, 
  typename PrecT
> 
class TrilinosLinearIterativeSolver 
  : public LinearIterativeSolverBase<
      TrilinosLinearIterativeSolver<
        SolverT,
        PrecT
      > 
    >
{

  private:

    friend LinearIterativeSolvers;
    typedef LinearIterativeSolverBase<TrilinosLinearIterativeSolver<SolverT, PrecT>> base_type;  
    typedef Epetra_RowMatrix matrix_type;


  public: 

    template <typename T>
    void resetLinearSystem(const core::Matrix<T>& A) {
      static_assert(std::is_base_of<Epetra_RowMatrix, std::decay_t<T>>::value, "Error: the supplied linear system cannot be used with the linear solver due to type incompatibility");
      solver->SetUserMatrix(A.data(), true);
    }


    template <typename T, typename U>
    void resetLinearSystem(const core::Matrix<T>& A, const core::Matrix<U>& P) {
      assert(A.data()->NumMyRows() == P.data()->NumMyRows());
      static_assert(std::is_base_of<Epetra_RowMatrix, std::decay_t<T>>::value, "Error: the supplied linear system cannot be used with the linear solver due to type incompatibility");
      static_assert(std::is_base_of<Epetra_RowMatrix, std::decay_t<U>>::value, "Error: the supplied preconditioner cannot be used with the linear solver due to type incompatibility");
      solver->SetUserMatrix(A.data());
      solver->SetPrecMatrix(P.data());
    }


    /**
     * @brief  Specify the matrix used to compute the preconditioner
     *
     * @param  P A preconditioner matrix 
     * @return void
     */
    template <typename T>
    void setPreconditionerMatrix(const core::Matrix<T>& P) {
      assert(solver->GetUserMatrix() != 0);
      assert(solver->GetUserMatrix()->NumMyRows() == P.data()->NumMyRows());
      static_assert(std::is_base_of<Epetra_RowMatrix, std::decay_t<T>>::value, "Error: the supplied preconditioner cannot be used with the linear solver due to type incompatibility");
      solver->SetPrecMatrix(P.data());
    }


    /**
     * @brief  Use the linear system matrix to compute the preconditioner
     *
     * @param  none
     * @return void
     */
    void setPreconditionerMatrix() {
      assert(solver->GetUserMatrix() != 0);
      solver->SetPrecMatrix(solver->GetUserMatrix()); 
    }


    template <typename T>
    auto solve(const T& b) {
      if (solver->GetUserMatrix() == 0) {std::cerr << "Error: the linear system has not been assigned yet" << std::endl; assert(0);}
      static_assert(core::details::vector_traits<T>::is_trilinos, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");

      T x(b.data()->Map());
      solver->SetLHS(&x);
      solver->SetRHS(b.data());
      solver->Iterate(maxIters_, tolerance_); 
      return x;
    }


    template <typename T, typename U>
    void solve(const T& b, U& x) {
      assert(solver->GetUserMatrix() != 0);
      static_assert(core::details::vector_traits<T>::is_trilinos, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");
      static_assert(core::details::vector_traits<U>::is_trilinos, "Error: the supplied result vector cannot be used with the linear solver due to type incompatibility");

      solver->SetLHS(&x);
      solver->SetRHS(&b);
      solver->Iterate(maxIters_, tolerance_);
    }


    int getMaxIterations() {
      return maxIters_;
    }


    void setMaxIterations(int maxIters) {
      maxIters_ = maxIters;
    }


    double getTolerance() {
      return tolerance_;
    }
 

    void setTolerance(double tolerance) {
      tolerance_ = tolerance;
    }


  private:

    TrilinosLinearIterativeSolver() = delete;
    TrilinosLinearIterativeSolver(const Epetra_RowMatrix& A) : 
      base_type(), solver(std::make_unique<AztecOO>()), maxIters_(0), tolerance_(1.0e-5)
    {
      resetLinearSystem(A);
      SolverT::set_solver_type(solver);
      PrecT::set_preconditioner_type(solver);
    }


  private:

    std::unique_ptr<AztecOO> solver;
 
    int maxIters_;
    double tolerance_; 
};

} //end namespace solvers

#endif
