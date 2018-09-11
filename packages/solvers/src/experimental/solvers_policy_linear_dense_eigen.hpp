
#ifndef SOLVERS_EXPERIMENTAL_POLICY_LINEAR_DENSE_EIGEN_HPP
#define SOLVERS_EXPERIMENTAL_POLICY_LINEAR_DENSE_EIGEN_HPP

#include <Eigen/Core>


namespace solvers {

// Policy to solve a dense linear system defined by an Eigen matrix
template <
  typename SolverT,
  typename MatrixT
>
struct SolversLinearDenseEigenPolicy {

  static void resetLinearSystem(std::shared_ptr<SolverT>& solver, const MatrixT& A) {
    solver->compute(*A.data());
  }


  template <typename VectorT>
  static auto solve(
    std::shared_ptr<SolverT> solver,
    const VectorT& b
  ) {
    return VectorT(solver->solve(*b.data()));
  }
};

} // end namespace solvers

#endif
