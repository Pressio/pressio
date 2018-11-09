#ifndef SOLVERS_LINEAR_WRAPPER_EIGEN_HPP
#define SOLVERS_LINEAR_WRAPPER_EIGEN_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../../../core/src/meta/core_meta_detection_idiom.hpp"


namespace rompp {
namespace solvers {


template <typename SolverT>
class SolversLinearDirectWrapperEigen {

  public:

    SolversLinearDirectWrapperEigen() : solver_() {}
    virtual ~SolversLinearDirectWrapperEigen() = default;


    /**
     * Update the matrix that describes the linear system.
     */
    template <
      typename MatrixT,
      typename core::meta::enable_if_t<
        core::details::traits<MatrixT>::wrapped_package_identifier == core::details::WrappedPackageIdentifier::Eigen
      >* = nullptr
    >
    void resetLinearSystem(const MatrixT& A) {
      solver_.compute(*A.data());
    }


    /**
     * Solve the linear system.
     */
    template <
      typename VectorT,
      typename core::meta::enable_if_t<
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
class SolversLinearIterativeWrapperEigen : public SolversLinearDirectWrapperEigen<SolverT> {


  public:

    SolversLinearIterativeWrapperEigen() : SolversLinearDirectWrapperEigen<SolverT>() {}
    ~SolversLinearIterativeWrapperEigen() = default;


    /**
     * Set the maximum number of iterations.
     */

    void setMaxIterations(int maxIters) {
      this->solver_.setMaxIterations(maxIters);
    }


    /**
     * Set the algorithm tolerance.
     */
    void setTolerance(double tolerance) {
      this->solver_.setTolerance(tolerance);
    }

};


} // end namespace solvers
} // end namespace rompp

#endif
