
#ifndef SOLVERS_LINEAR_ITERATIVE_EIGEN_HPP_
#define SOLVERS_LINEAR_ITERATIVE_EIGEN_HPP_

#include <memory>

#include "solvers_ConfigDefs.hpp"
#include "solvers_forward_declarations.hpp"

#include "vector/core_vector_traits_exp.hpp"
#include "solvers_linear_iterative_base.hpp"


namespace solvers {

// Forward declarations
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
    typedef LinearIterativeSolverBase<EigenLinearIterativeSolver<SolverT>> base_type;
    typedef typename SolverT::MatrixType matrix_type; 


  public: 

    EigenLinearIterativeSolver(EigenLinearIterativeSolver&& other) : solver(std::move(other.solver)), rows_(other.rows_) {}
    

    void resetLinearSystem(const core::matrix<matrix_type>& A) {
      rows_ = A.rows();
      solver->compute(*A.data());
    }


    template <typename T>
    auto solve(const T& b) {
      static_assert(core::details::vector_traits<T>::is_eigen, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");
      assert(core::details::vector_traits<T>::is_dynamic || rows_ == core::details::vector_traits<T>::rows);   

      T x = T(solver->solve(*b.data()).rhs());
      return x;
    }


    template <typename T, 
      typename U
    >
    void solve(const T& b, 
      U& x
    ) {
      static_assert(core::details::vector_traits<T>::is_eigen, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");
      static_assert(core::details::vector_traits<U>::is_eigen, "Error: the supplied result vector cannot be used with the linear solver due to type incompatibiliy");
      assert(core::details::vector_traits<T>::is_dynamic || rows_ == core::details::vector_traits<T>::rows);
      assert(core::details::vector_traits<U>::is_dynamic || rows_ == core::details::vector_traits<U>::rows);

      typename core::details::vector_traits<U>::wrapped_type tmp = solver->solve(*b.data()).rhs();
      x = U(tmp);
    }


    int getMaxIterations() {
      return solver.maxIterations();
    }


    void setMaxIterations(int maxIters) {
      solver->setMaxIterations(maxIters);
    }


    double getTolerance() {
      return solver->tolerance();
    }

 
    void setTolerance(double tolerance) {
      solver->setTolerance(tolerance);
    }


  private:

    EigenLinearIterativeSolver() = delete;
    EigenLinearIterativeSolver(const core::matrix<matrix_type>& A) : 
      base_type(), solver(std::make_unique<SolverT>()), rows_(0)
    {
      resetLinearSystem(A);
    }


  private:
 
    std::unique_ptr<SolverT> solver;
    int rows_;
};

} //end namespace solvers

#endif
