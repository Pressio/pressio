
#ifndef SOLVERS_EXPERIMENTAL_POLICY_LINEAR_ITERATIVE_EIGEN_HPP
#define SOLVERS_EXPERIMENTAL_POLICY_LINEAR_ITERATIVE_EIGEN_HPP

#include <Eigen/Core>
#include <memory>


namespace solvers {

// Policy to solve a linear system
template <
  typename SolverT, 
  typename MatrixT
>
struct SolversLinearIterativeEigenPolicy {

  static void resetLinearSystem(std::shared_ptr<SolverT>& solver, const MatrixT& A) {
    solver->compute(*A.data());
  }


  template <typename VectorT>
  static auto solve(
    std::shared_ptr<SolverT> solver, 
    const VectorT& b, 
    int maxIters, 
    double tolerance
  ) {
    solver->setMaxIterations(maxIters);
    solver->setTolerance(tolerance);
    return VectorT(solver->solve(*b.data()));
  }
};

} // end namespace solvers

#endif