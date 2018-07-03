
#ifndef SOLVERS_TRILINOS_LINEAR_ITERATIVE_HPP_
#define SOLVERS_TRILINOS_LINEAR_ITERATIVE_HPP_

#include "solvers_ConfigDefs.hpp"
#include "solvers_forward_declarations.hpp"

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
    typedef Epetra_RowMatrix matrix_type;


  public: 

    void resetLinearSystem(const matrix_type& A) {
      solver.SetUserMatrix(&A, true);
    }


    void resetLinearSystem(const matrix_type& A, const matrix_type& P) {
      assert(A.NumMyRows() == P.NumMyRows());
      solver.SetUserMatrix(&A);
      solver.SetPrecMatrix(&P);
    }


    /**
     * @brief  Specify the matrix used to compute the preconditioner
     *
     * @param  P A preconditioner matrix 
     * @return void
     */
    void setPreconditionerMatrix(const matrix_type& P) {
      assert(solver.GetUserMatrix() != 0);
      assert(solver.GetUserMatrix()->NumMyRows() == P.NumMyRows());
      solver.SetPrecMatrix(&P);
    }


    /**
     * @brief  Use the linear system matrix to compute the preconditioner
     *
     * @param  none
     * @return void
     */
    void setPreconditionerMatrix() {
      assert(solver.GetUserMatrix() != 0);
      solver.SetPrecMatrix(solver.GetUserMatrix()); 
    }


    template <typename T>
    auto solve(const T& b) {
      assert(solver.GetUserMatrix() != 0, "Error: the linear system has not been assigned yet");
      static_assert(core::vector_traits<T>::is_trilinos, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");

      T x(b.Map());
      solver.SetLHS(&x);
      solver.SetRHS(&b);
      solver.Iterate(maxIters_, tolerance_); 
      return x;
    }


    template <typename T, typename U>
    void solve(const T& b, U& x) {
      assert(solver.GetUserMatrix() != 0);
      static_assert(core::vector_traits<T>::is_trilinos, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");
      static_assert(core::vector_traits<U>::is_trilinos, "Error: the supplied result vector cannot be used with the linear solver due to type incompatibility");

      solver.SetLHS(&x);
      solver.SetRHS(&b);
      solver.Iterate(maxIters_, tolerance_);
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
      base()
    {
      resetLinearSystem(A);
      SolverT::set_solver_type(solver);
      PrecT::set_preconditioner_type(solver);
    }


  private:

    AztecOO solver;
 
    int maxIters_;
    double tolerance_; 
};

} //end namespace experimental
} //end namespace solvers

#endif
