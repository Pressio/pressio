
#ifdef HAVE_TRILINOS
#ifndef QR_TPETRA_MULTI_VECTOR_HPP_
#define QR_TPETRA_MULTI_VECTOR_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_forward_declarations.hpp"
#include "qr_solver_base.hpp"

#include "../../CORE_ALL"
#include <Eigen/OrderingMethods>
#include<Eigen/SparseQR>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#include "AnasaziTsqrOrthoManager.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziSolverUtils.hpp"
#include "AnasaziTpetraAdapter.hpp"


namespace rompp{ namespace qr{

// using anasazi TSQR
template<typename matrix_type,
	 template <typename...> class Q_type,
	 typename R_type>
class QRSolver<matrix_type,
	       Q_type,
	       R_type,
	       ::rompp::qr::TSQR,
	       typename
	       std::enable_if<
		 core::meta::is_tpetra_multi_vector_wrapper<matrix_type>::value and
		 core::meta::is_core_matrix_wrapper<R_type>::value and
		 core::details::traits<R_type>::is_shared_mem and
		 core::details::traits<R_type>::is_dense
		 >::type
	       >
  : public QRSolverBase<QRSolver<matrix_type, Q_type, R_type, ::rompp::qr::TSQR>,
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

  using this_t = QRSolver<matrix_type, Q_type, R_type, ::rompp::qr::TSQR>;
  using base_t = QRSolverBase<this_t, Q_t, R_type,  matrix_type>;

  using MVTraits = Anasazi::MultiVecTraits<sc_t, MV>;
  //  using OP = Tpetra::Operator<sc_t>;
  //  using ortho_t = Anasazi::TsqrMatOrthoManager<sc_t, MV, OP>;
  using ortho_t = Anasazi::TsqrOrthoManager<sc_t, MV>;

  using mat_type = Teuchos::SerialDenseMatrix<int, sc_t>;
  using mat_ptr = Teuchos::RCP<mat_type>;

  const std::string label_ = "Anasazi";

public:
  QRSolver() {
    OM_ = std::make_shared<ortho_t>(label_);
  }

  ~QRSolver() = default;

private:

  // non-const because Anasazi::TSQR can take a non-cost matrix
  // but here we do not modify it
  void computeImpl(matrix_type & A){
    // A is a rompp:core::MultiVector<MV>

    // get number of cols
    auto n = A.globalNumVectors();
    // the row map
    auto ArowMap = A.getRCPDataMap();

    // // clone the incoming matrix because TSQR modifies the matrix
    // Teuchos::RCP<MV> A1 = MVTraits::CloneCopy( *A.data() );
    // // Tpetra::MatrixMarket::Writer<MV>::writeDense(std::cout,
    // // 						 *A1, "void",
    // "dfdfd");

    // create Ortho manager if not already existing
    if (!OM_){
      OM_ = std::make_shared<ortho_t>(label_);
    }

    // create Q factor
    if (!Qmat_)
      Qmat_ = std::make_shared<Q_t>(ArowMap, n);

    // create B
    mat_ptr B = Teuchos::rcp(new mat_type(n,n) );

    const int initialX1Rank = OM_->normalizeOutOfPlace(*A.data(),
						       *Qmat_->data(),
						       B);

    // store R factor
    if (!Rmat_)
      Rmat_ = std::make_shared<R_type>(B->values());

    auto err = OM_->orthonormError(*Qmat_->data());
    std::cout << " Rank = "
	      << initialX1Rank << " "
	      <<  n << " "
	      << " error = " << err
	      << std::endl;

    //    B->print(std::cout);

  }//end method

  const Q_t & cRefQFactorImpl() const {
    return *Qmat_;
  }//end method

  const R_type & cRefRFactorImpl() const {
    return *Rmat_;
  }//end method

private:
  friend base_t;

  std::shared_ptr<Q_t> Qmat_     = nullptr;
  std::shared_ptr<R_type> Rmat_  = nullptr;
  std::shared_ptr< ortho_t > OM_ = nullptr;

};//end class









// overload for:
// the input data is a wrapper of an Tpetra_MultiVector
// the Q factor is a Q_type wrapper of an Tpetra_MultiVector
// the R factor is stored into whatever it is passed but the
// R_type has to be a core matrix wrapper sharedmem
template<typename matrix_type,
	 template <typename...> class Q_type,
	 typename R_type>
class QRSolver<matrix_type,
	       Q_type,
	       R_type,
	       ::rompp::qr::Hacked,
	       typename
	       std::enable_if<
		 core::meta::is_tpetra_multi_vector_wrapper<matrix_type>::value and
		 core::meta::is_core_matrix_wrapper<R_type>::value and
		 core::details::traits<R_type>::is_shared_mem
		 >::type
	       >{
private:
  using MV = typename core::details::traits<matrix_type>::wrapped_t;
  using sc_t = typename core::details::traits<matrix_type>::scalar_t;
  using LO_t = typename core::details::traits<matrix_type>::local_ordinal_t;
  using GO_t = typename core::details::traits<matrix_type>::global_ordinal_t;
  using map_t = typename core::details::traits<matrix_type>::data_map_t;
  using node_t = typename core::details::traits<matrix_type>::node_t;
  using hexsp = typename core::details::traits<matrix_type>::host_exec_space_t;
  using Q_t = Q_type<MV>;

public:
  QRSolver() = default;
  ~QRSolver() = default;

  void compute(const matrix_type & A){
    auto m = A.globalLength();
    auto n = A.globalNumVectors();
    auto ArowMap = A.getRCPDataMap();
    //auto comm = A.comm();
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
    auto Qm = eQR.householderQ() * Eigen::MatrixXd::Identity(m,m);
    auto & Rm = eQR.matrixQR().template triangularView<Eigen::Upper>();

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

    // store R factor
    Rmat_ = std::make_shared<R_type>(Rm);

  }//end method
  //-----------------------------------------------------

  const Q_t & cRefQFactor() const {
    return *Qmat_;
  }//end method
  //-----------------------------------------------------

  const R_type & cRefRFactor() const {
    return *Rmat_;
  }//end method
  //-----------------------------------------------------

private:
  std::shared_ptr<Q_t> Qmat_;
  std::shared_ptr<R_type> Rmat_;
};//end class


}} // end namespace rompp::qr
#endif
#endif
