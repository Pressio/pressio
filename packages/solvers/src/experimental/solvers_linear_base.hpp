
#ifndef SOLVERS_EXPERIMENTAL_LINEAR_BASE_HPP
#define SOLVERS_EXPERIMENTAL_LINEAR_BASE_HPP

#include <memory>
#include <type_traits>

#include "solvers_meta_static_checks.hpp"


namespace rompp{
namespace solvers{


/**
 * Base class for linear solver implemented through CRTP.
 *
 * @section DESCRIPTION
 *
 * This class defines the public interface for a linear solver class.
 * Objects of the class cannot be created directly. To create a solver,
 * use the factory class LinearSolvers.
 */
template<
  typename SolverT,
  typename MatrixT,
  typename Derived
>
class LinearSolverBase {

public:

  /**
   * Initialize a new linear solver.
   *
   * @param A matrix representing the linear system to be solved.
   */
  template <
  typename CompatibleMatrixT,
  typename std::enable_if<
    solvers::meta::are_matrix_compatible<
      MatrixT,
      CompatibleMatrixT
      >::value, MatrixT*
    >::type = nullptr
  >
  void resetLinearSystem(const CompatibleMatrixT& A) {
    solver_->resetLinearSystem(A);
  }


  /**
   * Solve the linear system
   *
   * @param b RHS vector
   * @return solution vector
   */
  template <
    typename VectorLT,
    typename std::enable_if<
      solvers::meta::are_vector_matrix_compatible<
	VectorLT,
	MatrixT
        >::value, MatrixT*
      >::type = nullptr
    >
  VectorLT solve(const VectorLT& b){
    return this->underlying()._solve(b);
  }


  /**
   * Specify and solve the linear system
   *
   * @param A matrix representing the linear system to be solved
   * @param b RHS vector
   * @return solution vector
   */
  template <
    typename CompatibleMatrixT,
    typename VectorRT,
    typename std::enable_if<
      solvers::meta::are_vector_matrix_compatible<
	VectorRT,
	CompatibleMatrixT
        >::value,
      CompatibleMatrixT*
      >::type = nullptr
    >
  auto solve(const CompatibleMatrixT& A, const VectorRT& b)
    -> decltype(this->solve(b)){
    this->resetLinearSystem(A);
    return this->solve(b);
  }


  /**
   * @brief  Solve the linear system
   *
   * @param  b is the RHS vector
   * @param  x is the solution vector
   * @return void
   */
  template <
    typename VectorLT,
    typename VectorRT,
    typename std::enable_if<
      solvers::meta::are_vector_compatible<
	VectorLT,
	VectorRT
        >::value,
      VectorLT*
      >::type = nullptr
    >
  void solve(const VectorLT& b, VectorRT& x) {
    x = VectorRT(*this->solve(b).data());
  }


  /**
   * @brief  Specify and solve the linear system
   *
   * @param  A is the system matrix
   * @param  b is the RHS vector
   * @param  x is the solution vector
   * @return void
   */
  template <
    typename CompatibleMatrixT,
    typename VectorLT,
    typename VectorRT
    >
  void solve(const CompatibleMatrixT& A, const VectorLT& b, VectorRT& x) {
    this->resetLinearSystem(A);
    this->solve(b, x);
  }


protected:

  LinearSolverBase() : solver_(nullptr) {};


  LinearSolverBase(std::shared_ptr<SolverT> solver) : solver_(solver) {}


  LinearSolverBase(LinearSolverBase&& other) : solver_(std::move(other.solver_)) {}


  LinearSolverBase(const LinearSolverBase&) = delete;


  virtual ~LinearSolverBase() = default;


  std::shared_ptr<SolverT> getSolver() {
    return solver_;
  }


private:

  Derived& underlying() {
    return static_cast<Derived&>(*this);
  }


  Derived const& underlying() const {
    return static_cast<Derived const&>(*this);
  }


private:

  std::shared_ptr<SolverT> solver_;
};

} //end namespace solvers

}//end namespace rompp
#endif
