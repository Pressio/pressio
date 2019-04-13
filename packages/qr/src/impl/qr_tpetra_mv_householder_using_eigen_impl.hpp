
#ifdef HAVE_TRILINOS
#ifndef QR_TPETRA_MV_HOUSEHOLDER_USING_EIGEN_IMPL_HPP_
#define QR_TPETRA_MV_HOUSEHOLDER_USING_EIGEN_IMPL_HPP_

#include "../qr_rfactor_solve_impl.hpp"

#include <Eigen/OrderingMethods>
#include <Eigen/SparseQR>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

namespace rompp{ namespace qr{ namespace impl{

template<typename matrix_t, typename R_t, int n, int m,
	 template <typename...> class Q_type>
class TpetraMVHouseholderUsingEigen{

  using MV = typename core::details::traits<matrix_t>::wrapped_t;
  using sc_t = typename core::details::traits<matrix_t>::scalar_t;
  using LO_t = typename core::details::traits<matrix_t>::local_ordinal_t;
  using GO_t = typename core::details::traits<matrix_t>::global_ordinal_t;
  using map_t = typename core::details::traits<matrix_t>::data_map_t;
  using node_t = typename core::details::traits<matrix_t>::node_t;
  using hexsp = typename core::details::traits<matrix_t>::host_exec_space_t;

  using Q_t = Q_type<MV>;
  using eig_dyn_mat	= Eigen::MatrixXd;
  using eig_mat_w	= core::Matrix<eig_dyn_mat>;
  using help_impl_t	= QRHouseholderDenseEigenMatrixWrapper<
				eig_mat_w, R_t, n, m, Q_type>;
  help_impl_t myImpl_	= {};

public:
  TpetraMVHouseholderUsingEigen() = default;
  ~TpetraMVHouseholderUsingEigen() = default;

  template < typename vector_in_t, typename vector_out_t>
  void project(const vector_in_t & vecIn,
  	       vector_out_t & vecOut) const{
    core::ops::dot( *this->Qmat_, vecIn, vecOut );
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
    auto ArowMap = A.getRCPDataMap();
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
      Teuchos::rcp (new Teuchos::MpiComm<int> (MPI_COMM_SELF));

    // convert it to replicated eptra matrix
    using local_map_t = Tpetra::Map<LO_t, GO_t, node_t>;
    using rcp_local_map_t = Teuchos::RCP<const local_map_t>;
    rcp_local_map_t rcp_local_map = Teuchos::rcp( new local_map_t(rows, 0, comm) );

    using import_t = Tpetra::Import<LO_t, GO_t, node_t>;
    import_t importer(ArowMap, rcp_local_map);
    matrix_t A2(rcp_local_map, cols);
    A2.data()->doImport(*A.data(), importer, Tpetra::INSERT);

    // store it into an Eigen matrix
    core::Matrix<Eigen::MatrixXd> eA2W(rows,cols);
    for (int j=0;j<cols;j++){
      auto colData = A2.data()->getData(j);
      for (int i=0;i<rows;i++)
    	eA2W(i,j) = colData[i];
    }

    myImpl_.computeThinOutOfPlace(eA2W);

    // // do QR in Eigen
    // Eigen::HouseholderQR<Eigen::MatrixXd> eQR(*eA2W.data());
    // auto Qm = eQR.householderQ() * Eigen::MatrixXd::Identity(m,n);
    // auto & Rm = eQR.matrixQR().template triangularView<Eigen::Upper>();

    // // store R factor
    // //auto RFn = Rm.block(0,0,n,n);
    // Rmat_ = std::make_shared<R_type>(Rm);

    // store Q into replicated Tpetra::Multivector
    const auto & Q2 = *(myImpl_.getCRefQFactor().data());
    Q_t locQ( rcp_local_map, Q2.cols() );
    auto trilD = locQ.data();
    trilD->template sync<Kokkos::HostSpace>();

    auto v2d = trilD->template getLocalView<Kokkos::HostSpace>();
    auto c0 = Kokkos::subview(v2d, Kokkos::ALL(), 0);
    // //we are going to change the host view
    trilD->template modify<Kokkos::HostSpace>();
    for (int i=0;i<Q2.rows();i++)
      for (int j=0;j<Q2.cols();j++)
    	v2d(i,j) = Q2(i,j);

    // import from local to distributed
    Qmat_ = std::make_shared<Q_t>(ArowMap, Q2.cols());
    import_t importer2(rcp_local_map, ArowMap);
    Qmat_->data()->doImport(*locQ.data(), importer2, Tpetra::INSERT);
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
