
#ifdef HAVE_TRILINOS
#ifndef QR_TPETRA_MULTI_VECTOR_HOUSEHOLDER_HPP_
#define QR_TPETRA_MULTI_VECTOR_HOUSEHOLDER_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_forward_declarations.hpp"
#include "qr_solver_base.hpp"

#include "../../CORE_ALL"
#include <Eigen/OrderingMethods>
#include<Eigen/SparseQR>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>


namespace rompp{ namespace qr{

// Tpetra multivector, householder
template<typename matrix_type,
	 typename R_type,
	 template <typename...> class Q_type>
class QRSolver<matrix_type,
	       R_type,
	       ::rompp::qr::Householder,
	       Q_type,
	       typename
	       std::enable_if<
		 core::meta::is_tpetra_multi_vector_wrapper<matrix_type>::value and
		 core::meta::is_core_matrix_wrapper<R_type>::value and
		 core::details::traits<R_type>::is_shared_mem and
		 core::details::traits<R_type>::is_dense
		 >::type
	       >
  : public QRSolverBase<QRSolver<matrix_type, R_type, ::rompp::qr::Householder, Q_type>,
			Q_type<typename core::details::traits<matrix_type>::wrapped_t>,
			R_type,
			matrix_type>{

  using MV = typename core::details::traits<matrix_type>::wrapped_t;
  using sc_t = typename core::details::traits<matrix_type>::scalar_t;
  using LO_t = typename core::details::traits<matrix_type>::local_ordinal_t;
  using GO_t = typename core::details::traits<matrix_type>::global_ordinal_t;
  using map_t = typename core::details::traits<matrix_type>::data_map_t;
  using node_t = typename core::details::traits<matrix_type>::node_t;
  using hexsp = typename core::details::traits<matrix_type>::host_exec_space_t;
  using Q_t = Q_type<MV>;

  using this_t = QRSolver<matrix_type, R_type, ::rompp::qr::Householder, Q_type>;
  using base_t = QRSolverBase<this_t, Q_t, R_type,  matrix_type>;

  using base_t::Qmat_;
  using base_t::Rmat_;

public:
  QRSolver() = default;
  ~QRSolver() = default;

private:

  void computeThinImpl(matrix_type & A){

    auto m = A.globalLength();
    auto n = A.globalNumVectors();
    auto ArowMap = A.getRCPDataMap();
    Teuchos::RCP<const Teuchos::Comm<int> > comm =
      Teuchos::rcp (new Teuchos::MpiComm<int> (MPI_COMM_SELF));

    // convert it to replicated eptra matrix
    using local_map_t = Tpetra::Map<LO_t, GO_t, node_t>;
    using rcp_local_map_t = Teuchos::RCP<const local_map_t>;
    rcp_local_map_t rcp_local_map = Teuchos::rcp( new local_map_t(m, 0, comm) );

    using import_t = Tpetra::Import<LO_t, GO_t, node_t>;
    import_t importer(ArowMap, rcp_local_map);
    matrix_type A2(rcp_local_map, n);
    A2.data()->doImport(*A.data(), importer, Tpetra::INSERT);

    // store it into an Eigen matrix
    core::Matrix<Eigen::MatrixXd> eA2W(m,n);
    for (int j=0;j<n;j++){
      auto colData = A2.data()->getData(j);
      for (int i=0;i<m;i++)
    	eA2W(i,j) = colData[i];
    }

    // do QR in Eigen
    Eigen::HouseholderQR<Eigen::MatrixXd> eQR(*eA2W.data());
    auto Qm = eQR.householderQ() * Eigen::MatrixXd::Identity(m,n);
    auto & Rm = eQR.matrixQR().template triangularView<Eigen::Upper>();

    // store R factor
    //auto RFn = Rm.block(0,0,n,n);
    Rmat_ = std::make_shared<R_type>(Rm);

    // store Q into replicated Tpetra::Multivector
    Q_t locQ( rcp_local_map, Qm.cols() );
    auto trilD = locQ.data();
    trilD->template sync<Kokkos::HostSpace>();

    auto v2d = trilD->template getLocalView<Kokkos::HostSpace>();
    auto c0 = Kokkos::subview(v2d, Kokkos::ALL(), 0);
    // //we are going to change the host view
    trilD->template modify<Kokkos::HostSpace>();
    for (int i=0;i<Qm.rows();i++)
      for (int j=0;j<Qm.cols();j++)
    	v2d(i,j) = Qm(i,j);

    // import from local to distributed
    Qmat_ = std::make_shared<Q_t>(ArowMap, Qm.cols());
    import_t importer2(rcp_local_map, ArowMap);
    Qmat_->data()->doImport(*locQ.data(), importer2, Tpetra::INSERT);

  }

  const Q_t & cRefQFactorImpl() const {
    return *Qmat_;
  }

  const R_type & cRefRFactorImpl() const {
    return *Rmat_;
  }

  friend base_t;

};//end class

}} // end namespace rompp::qr
#endif
#endif
