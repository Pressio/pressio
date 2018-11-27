
#ifndef QR_HACKED_HPP_
#define QR_HACKED_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_forward_declarations.hpp"
#include "../../CORE_ALL"
#include <Eigen/OrderingMethods>
#include<Eigen/SparseQR>
#ifdef HAVE_TRILINOS
#include <Epetra_Import.h>
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Teuchos_DefaultMpiComm.hpp>
#endif

namespace rompp{ namespace qr{ namespace hack{


// overload for:
// the input data is a wrapper for eigen DENSE matrix
// the Q factor is a Q_type wrapper of an eigen dense matrix
// the R factor is stored into whatever it is passed but the
// R_type has to be a core matrix wrapper sharedmem
template<typename matrix_type,
	 template <typename...> class Q_type,
	 typename R_type>
class QRSolver<matrix_type, Q_type, R_type,
	       core::meta::enable_if_t<
		 core::meta::is_eigen_dense_matrix_wrapper<matrix_type>::value and
		 core::meta::is_core_matrix_wrapper<R_type>::value and
		 core::details::traits<R_type>::is_shared_mem
		 >
	       >{
private:
  using sc_t = typename core::details::traits<matrix_type>::scalar_t;
  using Q_t = Q_type<Eigen::MatrixXd>;
public:
  QRSolver() = default;
  ~QRSolver() = default;

  void compute(const matrix_type & A){
    using native_mat_type = typename core::details::traits<matrix_type>::wrapped_t;
    /// do QR in Eigen
    Eigen::HouseholderQR<native_mat_type> eQR(*A.data());
    auto Qm = eQR.householderQ() * Eigen::MatrixXd::Identity(A.rows(), A.rows());
    auto & Rm = eQR.matrixQR().template triangularView<Eigen::Upper>();
    Qmat_ = std::make_shared<Q_t>(Qm);
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





// overload for:
// the input data is a wrapper for eigen SPARSE matrix
// the Q factor is a Q_type wrapper of an eigen dense matrix
// the R factor is stored into whatever it is passed but the
// R_type has to be a core matrix wrapper sharedmem
template<typename matrix_type,
	 template <typename...> class Q_type,
	 typename R_type>
class QRSolver<matrix_type, Q_type, R_type,
	       core::meta::enable_if_t<
		 core::meta::is_eigen_sparse_matrix_wrapper<matrix_type>::value and
		 core::meta::is_core_matrix_wrapper<R_type>::value and
		 core::details::traits<R_type>::is_shared_mem
		 >
	       >{
private:
  using sc_t = typename core::details::traits<matrix_type>::scalar_t;
  using Q_t = Q_type<Eigen::MatrixXd>;
public:
  QRSolver() = default;
  ~QRSolver() = default;

  void compute(const matrix_type & A){

    // I cannot get the spase QR to give me the Q matrix so I use for now dense QR

    //    using native_mat_type = typename core::details::traits<matrix_type>::wrapped_t;
    // using ord_type = typename core::details::traits<matrix_type>::ordinal_t;
    // assert(A.isCompressed());
    // using ordering_t = Eigen::NaturalOrdering<ord_type>;
    // Eigen::SparseQR<native_mat_type, ordering_t> eQR(*A.data());
    // assert(eQR.info() == Eigen::Success);
    // auto & Rm = eQR.matrixR().template triangularView<Eigen::Upper>();
    // auto GG = eQR.matrixQ();
    // Eigen::SparseMatrix<double, Eigen::RowMajor> sid;
    // GG.evalTo(sid);

    Eigen::MatrixXd Ad(*A.data());
    Eigen::HouseholderQR<Eigen::MatrixXd> eQR(Ad);
    auto Qm = eQR.householderQ() * Eigen::MatrixXd::Identity(A.rows(), A.rows());
    auto & Rm = eQR.matrixQR().template triangularView<Eigen::Upper>();
    Qmat_ = std::make_shared<Q_t>(Qm);
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





#ifdef HAVE_TRILINOS
// overload for:
// the input data is a wrapper of an Epetra_MultiVector
// the Q factor is a Q_type wrapper of an Epetra_MultiVector
// the R factor is stored into whatever it is passed but the
// R_type has to be a core matrix wrapper sharedmem
template<typename matrix_type,
	 template <typename...> class Q_type,
	 typename R_type>
class QRSolver<matrix_type, Q_type, R_type,
	       typename
	       std::enable_if<
     core::meta::is_core_multi_vector_wrapper<matrix_type>::value and
		 core::meta::is_multi_vector_epetra<
		   typename core::details::traits<matrix_type>::wrapped_t
		   >::value and
		 core::meta::is_core_matrix_wrapper<R_type>::value and
		 core::details::traits<R_type>::is_shared_mem
		 >::type
	       >{
private:
  using MV = Epetra_MultiVector;
  using sc_t = typename core::details::traits<matrix_type>::scalar_t;
  using Q_t = Q_type<MV>;

public:
  QRSolver() = default;
  ~QRSolver() = default;


  void compute(const matrix_type & A){
    auto m = A.globalLength();
    auto n = A.globalNumVectors();
    auto & ArowMap = A.getDataMap();

    // convert it to replicated eptra matrix
    Epetra_LocalMap locMap(m, 0, A.commCRef());
    Epetra_Import importer(locMap, ArowMap);
    matrix_type A2(locMap, n);
    A2.data()->Import(*A.data(), importer, Insert);

    // store it into an Eigen matrix
    core::Matrix<Eigen::MatrixXd> eA2W(m,n);
    for (int i=0;i<m;i++)
      for (int j=0;j<n;j++)
    	eA2W(i,j) = A2(i,j);

    //std::cout << *eA2W.data();
    // do QR in Eigen
    Eigen::HouseholderQR<Eigen::MatrixXd> eQR(*eA2W.data());
    auto Qm = eQR.householderQ() * Eigen::MatrixXd::Identity(m,m);
    auto & Rm = eQR.matrixQR().template triangularView<Eigen::Upper>();

    // store Q into replicated Epetra_Multivector
    Q_t locQ(locMap,Qm.cols());
    for (int i=0;i<Qm.rows();i++)
      for (int j=0;j<Qm.cols();j++)
    	locQ(i,j) = Qm(i,j);

    // import from local to distributed
    Qmat_ = std::make_shared<Q_t>(ArowMap, Qm.cols());
    Epetra_Import importer2(ArowMap, locMap);
    Qmat_->data()->Import(*locQ.data(), importer2, Insert);
    //    Qmat_->data()->Print(std::cout);

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





// overload for:
// the input data is a wrapper of an Tpetra_MultiVector
// the Q factor is a Q_type wrapper of an Tpetra_MultiVector
// the R factor is stored into whatever it is passed but the
// R_type has to be a core matrix wrapper sharedmem
template<typename matrix_type,
	 template <typename...> class Q_type,
	 typename R_type>
class QRSolver<matrix_type, Q_type, R_type,
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

#endif

}}} // end namespace rompp::qr::hack
#endif
