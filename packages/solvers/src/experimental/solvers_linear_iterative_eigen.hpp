
#ifndef SOLVERS_LINEAR_ITERATIVE_EIGEN_HPP_
#define SOLVERS_LINEAR_ITERATIVE_EIGEN_HPP_

#include "solvers_ConfigDefs.hpp"
#include "solvers_forward_declarations.hpp"

namespace solvers {

// Fwd declaration
struct LinearIterativeSolvers;
 

/**
 * @brief Class that implements LinearIterativeSolverBase for Eigen
 */ 
template<typename SolverT>  
class EigenLinearIterativeSolver 
  : public LinearIterativeSolverBase<
      EigenLinearIterativeSolver<
        SolverT
      >
    >
{

  private:

    friend LinearIterativeSolvers;  
    typedef SolverT::MatrixType matrix_type; 


  public: 

    void resetLinearSystem(const matrix_type& A) {
      rows_ = A.rows();
      solver.compute(A);
    }


    template <typename T>
    auto solve(const T& b) {
      static_assert(core::vector_traits<T>::is_eigen, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");
      assert(core::vector_traits<T>::is_dynamic || rows_ == core::vector_traits<T>::rows);   

      return solver.solve(b);
    }


    template <typename T, 
      typename U
    >
    void solve(const T& b, 
      U& x
    ) {
      static_assert(core::vector_traits<T>::is_eigen, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");
      static_assert(core::vector_traits<U>::is_eigen, "Error: the supplied result vector cannot be used with the linear solver due to type incompatibiliy");
      assert(core::vector_traits<T>::is_dynamic || rows_ == core::vector_traits<T>::rows);
      assert(core::vector<traits<U>::is_dynamic || rows_ == core::vector_traits<U>::rows);

      x = solver.solve(b)
    }


    int getMaxIterations() {
      return solver.maxIterations();
    }


    void setMaxIterations(int maxIters) {
      solver.setMaxIterations(maxIters);
    }


    double getTolerance() {
      return solver.tolerance();
    }

 
    void setTolerance(double tolerance) {
      solver.setTolerance(tolerance);
    }


  private:

    EigenLinearIterativeSolver() = delete;
    EigenLinearIterativeSolver(const MatrixType& A) : 
      base(), rows_(0)
    {
      resetLinearSystem(A);
    }


  private:
 
    SolvType solver;
    int rows_;
};

} //end namespace solvers

#endif
