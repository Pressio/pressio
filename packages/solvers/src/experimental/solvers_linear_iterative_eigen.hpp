
#ifndef SOLVERS_LINEAR_ITERATIVE_EIGEN_HPP_
#define SOLVERS_LINEAR_ITERATIVE_EIGEN_HPP_

<<<<<<< HEAD
#include <memory>

#include "solvers_ConfigDefs.hpp"
#include "solvers_forward_declarations.hpp"

#include "vector/core_vector_traits_exp.hpp"
#include "solvers_linear_iterative_base.hpp"


namespace solvers {

// Forward declarations
struct LinearIterativeSolvers;

=======
#include "solvers_ConfigDefs.hpp"
#include "solvers_forward_declarations.hpp"

namespace solvers {

// Fwd declaration
struct LinearIterativeSolvers;
 
>>>>>>> 8d54bd03f342d94a9ed90f26faae390adfe2ce0c

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

<<<<<<< HEAD
    friend LinearIterativeSolvers;
    typedef LinearIterativeSolverBase<EigenLinearIterativeSolver<SolverT>> base_type;
    typedef typename SolverT::MatrixType matrix_type; 
=======
    friend LinearIterativeSolvers;  
    typedef SolverT::MatrixType matrix_type; 
>>>>>>> 8d54bd03f342d94a9ed90f26faae390adfe2ce0c


  public: 

<<<<<<< HEAD
    EigenLinearIterativeSolver(EigenLinearIterativeSolver&& other) : solver(std::move(other.solver)), rows_(other.rows_) {}
    

    void resetLinearSystem(const core::matrix<matrix_type>& A) {
      rows_ = A.rows();
      solver->compute(*A.data());
=======
    void resetLinearSystem(const matrix_type& A) {
      rows_ = A.rows();
      solver.compute(A);
>>>>>>> 8d54bd03f342d94a9ed90f26faae390adfe2ce0c
    }


    template <typename T>
    auto solve(const T& b) {
<<<<<<< HEAD
      static_assert(core::details::vector_traits<T>::is_eigen, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");
      assert(core::details::vector_traits<T>::is_dynamic || rows_ == core::details::vector_traits<T>::rows);   

      T x = solver->solve(*b.data()).rhs();
      return x;
=======
      static_assert(core::vector_traits<T>::is_eigen, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");
      assert(core::vector_traits<T>::is_dynamic || rows_ == core::vector_traits<T>::rows);   

      return solver.solve(b);
>>>>>>> 8d54bd03f342d94a9ed90f26faae390adfe2ce0c
    }


    template <typename T, 
      typename U
    >
    void solve(const T& b, 
      U& x
    ) {
<<<<<<< HEAD
      static_assert(core::details::vector_traits<T>::is_eigen, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");
      static_assert(core::details::vector_traits<U>::is_eigen, "Error: the supplied result vector cannot be used with the linear solver due to type incompatibiliy");
      assert(core::details::vector_traits<T>::is_dynamic || rows_ == core::details::vector_traits<T>::rows);
      assert(core::details::vector_traits<U>::is_dynamic || rows_ == core::details::vector_traits<U>::rows);

      x = solver->solve(*b.data()).rhs();
=======
      static_assert(core::vector_traits<T>::is_eigen, "Error: the supplied RHS vector cannot be used with the linear solver due to type incompatibility");
      static_assert(core::vector_traits<U>::is_eigen, "Error: the supplied result vector cannot be used with the linear solver due to type incompatibiliy");
      assert(core::vector_traits<T>::is_dynamic || rows_ == core::vector_traits<T>::rows);
      assert(core::vector<traits<U>::is_dynamic || rows_ == core::vector_traits<U>::rows);

      x = solver.solve(b)
>>>>>>> 8d54bd03f342d94a9ed90f26faae390adfe2ce0c
    }


    int getMaxIterations() {
      return solver.maxIterations();
    }


    void setMaxIterations(int maxIters) {
<<<<<<< HEAD
      solver->setMaxIterations(maxIters);
=======
      solver.setMaxIterations(maxIters);
>>>>>>> 8d54bd03f342d94a9ed90f26faae390adfe2ce0c
    }


    double getTolerance() {
<<<<<<< HEAD
      return solver->tolerance();
=======
      return solver.tolerance();
>>>>>>> 8d54bd03f342d94a9ed90f26faae390adfe2ce0c
    }

 
    void setTolerance(double tolerance) {
<<<<<<< HEAD
      solver->setTolerance(tolerance);
=======
      solver.setTolerance(tolerance);
>>>>>>> 8d54bd03f342d94a9ed90f26faae390adfe2ce0c
    }


  private:

    EigenLinearIterativeSolver() = delete;
<<<<<<< HEAD
    EigenLinearIterativeSolver(const core::matrix<matrix_type>& A) : 
      base_type(), solver(std::make_unique<SolverT>()), rows_(0)
=======
    EigenLinearIterativeSolver(const MatrixType& A) : 
      base(), rows_(0)
>>>>>>> 8d54bd03f342d94a9ed90f26faae390adfe2ce0c
    {
      resetLinearSystem(A);
    }


  private:
 
<<<<<<< HEAD
    std::unique_ptr<SolverT> solver;
=======
    SolvType solver;
>>>>>>> 8d54bd03f342d94a9ed90f26faae390adfe2ce0c
    int rows_;
};

} //end namespace solvers

#endif
