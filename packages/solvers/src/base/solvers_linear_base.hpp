
#ifndef SOLVERS_LINEAR_BASE_HPP
#define SOLVERS_LINEAR_BASE_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../solvers_meta_static_checks.hpp"


namespace rompp{ namespace solvers{

/**
 * Base class for linear solver 
 *
 * @section DESCRIPTION
 *
 * This class defines the public interface for a linear solver class.
 * Objects of the class cannot be created directly. To create a solver,
 * use the factory class LinearSolvers.
 */
template<typename SolverT, typename MatrixT, typename Derived>
struct LinearBase {

  template <typename CompatibleMatrixT,
            typename std::enable_if<
              solvers::meta::are_matrix_compatible<
                MatrixT,
                CompatibleMatrixT
                >::value, MatrixT*
              >::type = nullptr
  >
  void resetLinearSystem(const CompatibleMatrixT& A) {
    this->underlying().resetLinearSystemImpl(A);
  }

  template <typename VectorT>
  void solve(const VectorT & b, VectorT& x) {
    this->underlying().solveImpl(b, x);
  }

  template <typename VectorT>
  void solve(const MatrixT & A, const VectorT & b, VectorT& x) {
    this->resetLinearSystem(A);
    this->solve(b, x);
  }

protected:
  LinearBase() = default;
  // LinearBase(std::shared_ptr<SolverT> solver) : solver_(solver) {}
  // LinearBase(LinearBase&& other) : solver_(std::move(other.solver_)) {}

  LinearBase(const LinearBase&) = delete;

  virtual ~LinearBase() = default;

  // std::shared_ptr<SolverT> getSolver() {
  //   return solver_;
  // }

private:
  // friend Derived;
  Derived& underlying(){ return static_cast<Derived&>(*this);  }
  Derived const& underlying() const { return static_cast<Derived const&>(*this); }

// private:
//   std::shared_ptr<SolverT> solver_ = nullptr;
};

}}//end namespace rompp::solvers
#endif









  // /**
  //  * @brief  Specify and solve the linear system
  //  *
  //  * @param  A is the system matrix
  //  * @param  b is the RHS vector
  //  * @param  x is the solution vector
  //  * @return void
  //  */
  // template <
  //   typename CompatibleMatrixT,
  //   typename VectorLT,
  //   typename VectorRT
  //   >
  // void solve(const CompatibleMatrixT& A, 
  //            const VectorLT& b, 
  //            VectorRT& x) {
  //   this->resetLinearSystem(A);
  //   // this->solve(b, x);
  //   this->underlying()._solve(b, x);
  // }

 //  /**
 //   * Solve the linear system
 //   *
 //   * @param b RHS vector
 //   * @return solution vector
 //   */
 //  template <
 //    typename VectorLT,
 //    typename std::enable_if<
 //      solvers::meta::are_vector_matrix_compatible<
  // VectorLT,
  // MatrixT
 //        >::value, MatrixT*
 //      >::type = nullptr
 //    >
 //  VectorLT solve(const VectorLT& b){
 //    return this->underlying()._solve(b);
 //  }


 //  *
 //   * Specify and solve the linear system
 //   *
 //   * @param A matrix representing the linear system to be solved
 //   * @param b RHS vector
 //   * @return solution vector
   
 //  template <
 //    typename CompatibleMatrixT,
 //    typename VectorRT,
 //    typename std::enable_if<
 //      solvers::meta::are_vector_matrix_compatible<
  // VectorRT,
  // CompatibleMatrixT
 //        >::value,
 //      CompatibleMatrixT*
 //      >::type = nullptr
 //    >
 //  auto solve(const CompatibleMatrixT& A, const VectorRT& b)
 //    -> decltype(this->solve(b)){
 //    this->resetLinearSystem(A);
 //    return this->solve(b);
 //  }
