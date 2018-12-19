
#ifdef HAVE_TRILINOS
#ifndef QR_EPETRA_MULTI_VECTOR_HOUSEHOLDER_HPP_
#define QR_EPETRA_MULTI_VECTOR_HOUSEHOLDER_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_forward_declarations.hpp"
#include "qr_solver_base.hpp"

#include "../../CORE_ALL"
#include <Eigen/OrderingMethods>
#include<Eigen/SparseQR>
#include <Epetra_Import.h>


namespace rompp{ namespace qr{

template<typename matrix_type,
	 typename R_type,
 	 template <typename...> class Q_type>
class QRSolver<matrix_type,
	       R_type,
	       ::rompp::qr::Householder,
	       Q_type,
	       typename
	       std::enable_if<
		 core::meta::is_epetra_multi_vector_wrapper<matrix_type>::value and
		 core::meta::is_core_matrix_wrapper<R_type>::value and
		 core::details::traits<R_type>::is_shared_mem and
		 core::details::traits<R_type>::is_dense
		 >::type
	       >
  : public QRSolverBase<QRSolver<matrix_type, R_type, ::rompp::qr::Householder, Q_type>,
			Q_type<typename core::details::traits<matrix_type>::wrapped_t>,
			R_type,
			matrix_type>{

private:
  using MV = Epetra_MultiVector;
  using sc_t = typename core::details::traits<matrix_type>::scalar_t;
  using Q_t = Q_type<MV>;

  using this_t = QRSolver<matrix_type, R_type, ::rompp::qr::Householder, Q_type>;
  using base_t = QRSolverBase<this_t, Q_t, R_type,  matrix_type>;
  friend base_t;

  using base_t::Qmat_;
  using base_t::Rmat_;

public:
  QRSolver() = default;
  ~QRSolver() = default;

private:

  void computeThinImpl(matrix_type & A){
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

    // do QR in Eigen
    Eigen::HouseholderQR<Eigen::MatrixXd> eQR(*eA2W.data());
    auto Qm = eQR.householderQ() * Eigen::MatrixXd::Identity(m,n);
    Eigen::MatrixXd Rm(eQR.matrixQR().template triangularView<Eigen::Upper>());

    // store R factor
    Rmat_ = std::make_shared<R_type>( Rm.block(0,0,n,n) );

    // store Q into replicated Epetra_Multivector
    Q_t locQ(locMap,Qm.cols());
    for (int i=0;i<Qm.rows();i++)
      for (int j=0;j<Qm.cols();j++)
    	locQ(i,j) = Qm(i,j);

    // import from local to distributed
    Qmat_ = std::make_shared<Q_t>(ArowMap, Qm.cols());
    Epetra_Import importer2(ArowMap, locMap);
    Qmat_->data()->Import(*locQ.data(), importer2, Insert);

  }

  const Q_t & cRefQFactorImpl() const {
    return *Qmat_;
  }

  const R_type & cRefRFactorImpl() const {
    return *Rmat_;
  }

};//end class


}} // end namespace rompp::qr
#endif
#endif
