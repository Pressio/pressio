#ifndef SOLVERS_LINEAR_WRAPPER_ARMADILLO_HPP
#define SOLVERS_LINEAR_WRAPPER_ARMADILLO_HPP

#include "../solvers_ConfigDefs.hpp"
#include "../../../containers/src/meta/containers_meta_detection_idiom.hpp"


namespace rompp {
namespace solvers {

template <typename MatrixT>
class SolversLinearDirectWrapperArmadillo {

  public:

    /**
     * Update the matrix that describes the linear system.
     */
    void resetLinearSystem(const MatrixT& other) {
      A = other;
    }


    /**
     * Solve the linear system
     */
    template <typename VectorT>
    VectorT solve(const VectorT& b) {
      return this->template _solveImpl<MatrixT>(b);
    }


  private:

    /**
     * Implement the solve method for a dense matrix.
     */
    template <
      typename DMatrixT,
      typename VectorT,
      typename ::rompp::mpl::enable_if_t<
        containers::details::traits<DMatrixT>::is_dense
      >* = nullptr
    >
    VectorT _solveImpl(const VectorT& b) const {
      return b;//VectorT(arma::solve(*A.data(), *b.data()));
    }


    /**
     * Implement the solve method for a sparse matrix.
     */
    template <
      typename DMatrixT,
      typename VectorT,
      typename ::rompp::mpl::enable_if_t<
        containers::details::traits<DMatrixT>::is_sparse
      >* = nullptr
    >
    VectorT _solveImpl(const VectorT& b) const {
      return b; // VectorT(arma::spsolve(*A.data(), *b.data()));
    }


  private:
    MatrixT A;
};


} // end namespace solvers
} // end namespace rompp

#endif
