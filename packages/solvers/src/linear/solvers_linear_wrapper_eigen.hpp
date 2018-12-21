#ifndef SOLVERS_LINEAR_WRAPPER_EIGEN_HPP
#define SOLVERS_LINEAR_WRAPPER_EIGEN_HPP

#include "../solvers_ConfigDefs.hpp"
#include <Eigen/Core>

namespace rompp { namespace solvers {


template <typename SolverT>
struct SolversLinearDirectWrapperEigen {

  SolversLinearDirectWrapperEigen() : solver_() {}
  virtual ~SolversLinearDirectWrapperEigen() = default;

  template <
    typename MatrixT,
    core::meta::enable_if_t<
      core::details::traits<MatrixT>::wrapped_package_identifier == core::details::WrappedPackageIdentifier::Eigen
    >* = nullptr
  >
  void resetLinearSystem(const MatrixT& A) {
    solver_.compute(*A.data());
  }

  template <
    typename VectorT,
    core::meta::enable_if_t<
      core::details::traits<VectorT>::wrapped_package_identifier == core::details::WrappedPackageIdentifier::Eigen
    >* = nullptr
  >
  VectorT solve(const VectorT& b) {
    return VectorT(solver_.solve(*b.data()));
  }

  protected:
    SolverT solver_;
};


template <typename SolverT>
struct SolversLinearIterativeWrapperEigen : public SolversLinearDirectWrapperEigen<SolverT> {

  SolversLinearIterativeWrapperEigen() : SolversLinearDirectWrapperEigen<SolverT>() {}
  ~SolversLinearIterativeWrapperEigen() = default;

  void setMaxIterations(int maxIters) {
    this->solver_.setMaxIterations(maxIters);
  }

  void setTolerance(double tolerance) {
    this->solver_.setTolerance(tolerance);
  }

};


}} // end namespace rompp::solvers

#endif
