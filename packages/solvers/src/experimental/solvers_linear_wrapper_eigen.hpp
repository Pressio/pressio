#ifndef SOLVERS_LINEAR_WRAPPER_EIGEN_HPP
#define SOLVERS_LINEAR_WRAPPER_EIGEN_HPP

#include "../solvers_ConfigDefs.hpp"

namespace rompp { namespace solvers {


template <typename SolverT>
struct SolversLinearDirectWrapperEigen {

  SolversLinearDirectWrapperEigen() : solver_() {}
  virtual ~SolversLinearDirectWrapperEigen() = default;

  template <
    typename MatrixT,
    typename ::rompp::mpl::enable_if_t<
      containers::details::traits<MatrixT>::wrapped_package_identifier == containers::details::WrappedPackageIdentifier::Eigen
    >* = nullptr
  >
  void resetLinearSystem(const MatrixT& A) {
    solver_.compute(*A.data());
  }

  template <
    typename VectorT,
    typename ::rompp::mpl::enable_if_t<
      containers::details::traits<VectorT>::wrapped_package_identifier == containers::details::WrappedPackageIdentifier::Eigen
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
