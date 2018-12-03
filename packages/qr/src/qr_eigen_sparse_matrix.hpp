
#ifndef QR_EIGEN_SPARSE_MATRIX_HPP_
#define QR_EIGEN_SPARSE_MATRIX_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_forward_declarations.hpp"
#include "../../CORE_ALL"
#include <Eigen/OrderingMethods>
#include<Eigen/SparseQR>

namespace rompp{ namespace qr{

// overload for:
// the input data is a wrapper for eigen SPARSE matrix
// the Q factor is a Q_type wrapper of an eigen dense matrix
// the R factor is stored into whatever it is passed but the
// R_type has to be a core matrix wrapper sharedmem
template<typename matrix_type,
	 template <typename...> class Q_type,
	 typename R_type>
class QRSolver<matrix_type,
	       Q_type,
	       R_type,
	       ::rompp::qr::Hacked,
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

    // using native_mat_type = typename core::details::traits<matrix_type>::wrapped_t;
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


}} // end namespace rompp::qr
#endif
