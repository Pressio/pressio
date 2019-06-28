
#ifdef HAVE_TRILINOS
#ifndef QR_EPETRA_MV_HOUSEHOLDER_USING_EIGEN_IMPL_HPP_
#define QR_EPETRA_MV_HOUSEHOLDER_USING_EIGEN_IMPL_HPP_

#include "../qr_rfactor_solve_impl.hpp"

#include <Eigen/OrderingMethods>
#include <Eigen/SparseQR>
#include <Epetra_Import.h>


namespace rompp{ namespace qr{ namespace impl{

template<typename matrix_t, typename R_t, int n, int m,
	 template <typename...> class Q_type>
class EpetraMVHouseholderUsingEigen{

  using MV = Epetra_MultiVector;
  using sc_t = typename algebra::details::traits<matrix_t>::scalar_t;
  using Q_t = Q_type<MV>;

  using eig_dyn_mat	= Eigen::MatrixXd;
  using eig_mat_w	= algebra::Matrix<eig_dyn_mat>;
  using help_impl_t	= QRHouseholderDenseEigenMatrixWrapper<
				eig_mat_w, R_t, n, m, Q_type>;
  help_impl_t myImpl_	= {};

public:
  EpetraMVHouseholderUsingEigen() = default;
  ~EpetraMVHouseholderUsingEigen() = default;

  template < typename vector_in_t, typename vector_out_t>
  void project(const vector_in_t & vecIn,
  	       vector_out_t & vecOut) const{
    algebra::ops::dot( *this->Qmat_, vecIn, vecOut );
  }

  template <typename vector_t>
  void doLinSolve(const vector_t & rhs, vector_t & y)const{
    myImpl_.template doLinSolve<vector_t>(rhs, y);
  }

  const Q_t & getCRefQFactor() const {
    return *this->Qmat_;
  }

  void computeThinOutOfPlace(matrix_t & A){
    auto rows = A.globalLength();
    auto cols = A.globalNumVectors();
    auto & ArowMap = A.getDataMap();

    // convert it to replicated eptra matrix
    Epetra_LocalMap locMap(rows, 0, A.commCRef());
    Epetra_Import importer(locMap, ArowMap);
    matrix_t A2(locMap, cols);
    A2.data()->Import(*A.data(), importer, Insert);

    // store it into an Eigen matrix
    eig_mat_w eA2W(rows,cols);
    for (int i=0;i<rows;i++)
      for (int j=0;j<cols;j++)
    	eA2W(i,j) = A2(i,j);

    myImpl_.computeThinOutOfPlace(eA2W);

    // // do QR in Eigen
    // hhold_t eQR(*eA2W.data());
    // auto Qm = eQR.householderQ() * eig_dyn_mat::Identity(rows,n);
    // eig_dyn_mat Rm(eQR.matrixQR().template triangularView<Eigen::Upper>());

    // // store R factor
    // Rmat_ = std::make_shared<R_type>( Rm.block(0,0,n,n) );

    // store Q into replicated Epetra_Multivector
    const auto & Q2 = *(myImpl_.getCRefQFactor().data());
    Q_t locQ(locMap,Q2.cols());
    for (int i=0;i<Q2.rows();i++)
      for (int j=0;j<Q2.cols();j++)
    	locQ(i,j) = Q2(i,j);

    // import from local to distributed
    Qmat_ = std::make_shared<Q_t>(ArowMap, Q2.cols());
    Epetra_Import importer2(ArowMap, locMap);
    Qmat_->data()->Import(*locQ.data(), importer2, Insert);
  }

private:

  // template <typename T = R_type,
  //       ::rompp::mpl::enable_if_t<
  //         !std::is_void<R_type>::value
  //         > * = nullptr>
  // const R_type & cRefRFactorImpl() const {
  //   return *Rmat_;
  // }

  // todo: these must be moved somewhere else
  mutable std::shared_ptr<Q_t> Qmat_	= nullptr;
  mutable std::shared_ptr<R_t> Rmat_	= nullptr;

};//end class


}}} // end namespace rompp::qr::impl
#endif
#endif
