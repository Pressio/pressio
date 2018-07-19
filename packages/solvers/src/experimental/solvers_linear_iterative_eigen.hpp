
#ifndef SOLVERS_LINEAR_ITERATIVE_EIGEN_HPP_
#define SOLVERS_LINEAR_ITERATIVE_EIGEN_HPP_

#include <memory>

#include "solvers_linear_iterative_base.hpp"
#include "matrix/core_matrix_traits_exp.hpp"
#include "vector/core_vector_traits_exp.hpp"


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
    

    /// Type of new linear system matrix is the same one used to initialize the solver  
    void resetLinearSystem(const core::Matrix<matrix_type>& A) {
      rows_ = A.rows();
      solver->compute(*A.data());
    }


    /// Scalar type of new linear system matrix differs from the one used to initialize the solver
    template <
      typename T,
      typename std::enable_if<
        !std::is_same<
          typename matrix_type::Scalar,
          typename core::details::matrix_traits<T>::wrapped_type::Scalar
        >::value, T
      >::type* = nullptr
    >
    void resetLinearSystem(const T& A) {
      static_assert(core::details::matrix_traits<T>::matrix_class == core::details::WrappedClass::Eigen, "Error: the supplied linear system matrix is not compatible with the linear solver"); 
      rows_ = A.rows();
      solver->compute(A.data()->template cast<typename matrix_type::Scalar>());
    }


    /// Scalar types of linear system matrix A and RHS vector b are the same
    template <
      typename T,
      typename std::enable_if<
        std::is_same<
          typename core::details::vector_traits<T>::wrapped_type::Scalar,
          typename matrix_type::Scalar
        >::value, T
      >::type* = nullptr
    >
    auto solve(const T& b) {
      typedef typename core::details::vector_traits<T> vector_traits;

      static_assert(vector_traits::vector_class == core::details::WrappedClass::Eigen, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");
      assert(vector_traits::is_dynamic || rows_ == vector_traits::rows);

      auto x = T(solver->solve(*b.data()));
      return x;
    }


    /// Scalar types of linear system matrix A and RHS vector b differs 
    template <
      typename T,
      typename std::enable_if<
        !std::is_same<
          typename core::details::vector_traits<T>::wrapped_type::Scalar, 
          typename matrix_type::Scalar
        >::value, T
      >::type* = nullptr 
    >
    auto solve(const T& b) {
      typedef typename core::details::vector_traits<T> vector_traits;

      static_assert(vector_traits::vector_class == core::details::WrappedClass::Eigen, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");
      assert(vector_traits::is_dynamic || rows_ == vector_traits::rows);   

      typedef typename vector_traits::wrapped_type::Scalar scalar_type;
      return T(solver->solve(b.data()->template cast<typename matrix_type::Scalar>()).template cast<scalar_type>()); 
    }


    /// Scalar types of RHS vector b and solution vector x are the smae
    template <typename T,
      typename U,
      typename std::enable_if<
        std::is_same<
          typename core::details::vector_traits<T>::wrapped_type::Scalar,
          typename core::details::vector_traits<U>::wrapped_type::Scalar
        >::value, T
      >::type* = nullptr
    >
    void solve(const T& b,
      U& x
    ) {
      typedef typename core::details::vector_traits<T> t_vector_traits;
      typedef typename core::details::vector_traits<U> u_vector_traits;

      static_assert(t_vector_traits::vector_class == core::details::WrappedClass::Eigen, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");
      static_assert(u_vector_traits::vector_class == core::details::WrappedClass::Eigen, "Error: the supplied result vector cannot be used with the linear solver due to type incompatibiliy");
      assert(t_vector_traits::is_dynamic || rows_ == t_vector_traits::rows);
      assert(u_vector_traits::is_dynamic || rows_ == u_vector_traits::rows);

      x= this->solve(b);
    }


    /// Scalar types of RHS vector b and solution vector x differs
    template <typename T, 
      typename U,
      typename std::enable_if<
        !std::is_same<
          typename core::details::vector_traits<T>::wrapped_type::Scalar,
          typename core::details::vector_traits<U>::wrapped_type::Scalar
        >::value, T
      >::type* = nullptr
    >
    void solve(const T& b, 
      U& x
    ) {
      typedef typename core::details::vector_traits<T> t_vector_traits;
      typedef typename core::details::vector_traits<U> u_vector_traits;

      static_assert(t_vector_traits::vector_class == core::details::WrappedClass::Eigen, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");
      static_assert(u_vector_traits::vector_class == core::details::WrappedClass::Eigen, "Error: the supplied result vector cannot be used with the linear solver due to type incompatibiliy");
      assert(t_vector_traits::is_dynamic || rows_ == t_vector_traits::rows);
      assert(u_vector_traits::is_dynamic || rows_ == u_vector_traits::rows);

      typedef typename u_vector_traits::wrapped_type::Scalar scalar_type;
      x = U(this->solve(b).data()->template cast<scalar_type>()); 
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
    EigenLinearIterativeSolver(const core::Matrix<matrix_type>& A) : 
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
