
#ifndef QR_EIGEN_SPARSE_MATRIX_HPP_
#define QR_EIGEN_SPARSE_MATRIX_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_forward_declarations.hpp"
#include "qr_solver_base.hpp"
#include "qr_rfactor_solve_impl.hpp"

#include "../../CORE_ALL"
#include <Eigen/OrderingMethods>
#include<Eigen/SparseQR>


namespace rompp{ namespace qr{

// the input data is a wrapper for eigen SPARSE matrix
template<typename matrix_type,
	 typename R_type,
	 template <typename...> class Q_type>
class QRSolver<matrix_type,
	       ::rompp::qr::Householder,
	       R_type,
	       Q_type,
	       core::meta::enable_if_t<
		 core::meta::is_eigen_sparse_matrix_wrapper<matrix_type>::value and
		 core::meta::is_core_matrix_wrapper<R_type>::value and
		 core::details::traits<R_type>::is_shared_mem and
		 core::details::traits<R_type>::is_dense
		 >
	       >
  : public QRSolverBase<QRSolver<matrix_type, ::rompp::qr::Householder, R_type, Q_type>,
			R_type, Q_type<Eigen::MatrixXd>, matrix_type>{

  using sc_t = typename core::details::traits<matrix_type>::scalar_t;
  using Q_t = Q_type<Eigen::MatrixXd>;

  using this_t = QRSolver<matrix_type, ::rompp::qr::Householder, R_type, Q_type>;
  using base_t = QRSolverBase<this_t, R_type, Q_t, matrix_type>;
  friend base_t;
  using base_t::Qmat_;
  using base_t::Rmat_;

public:
  QRSolver() = default;
  ~QRSolver() = default;

private:

  template <typename vector_in_t,
	    typename vector_out_t,
	    core::meta::enable_if_t<
	      core::meta::is_core_vector_wrapper<vector_in_t>::value and
	      core::meta::is_core_vector_wrapper<vector_out_t>::value and
	      // the type vector in should be from same package as Q
	      core::details::traits<vector_in_t>::wrapped_package_identifier ==
		core::details::traits<Q_t>::wrapped_package_identifier
	      > * = nullptr
	    >
  void projectImpl(const vector_in_t & vecIn,
		   vector_out_t & vecOut) const{
    core::ops::dot( *Qmat_, vecIn, vecOut );
  }


  void computeThinImpl(matrix_type & A){
    auto m = A.rows();
    auto n = A.cols();

    using native_mat_type = typename core::details::traits<matrix_type>::wrapped_t;
    using ord_type = typename core::details::traits<matrix_type>::ordinal_t;
    if (!A.isCompressed())
      A.compress();

    using ordering_t = Eigen::NaturalOrdering<ord_type>;
    Eigen::SparseQR<native_mat_type, ordering_t> eQR(*A.data());
    assert(eQR.info() == Eigen::Success);

    Qmat_ = std::make_shared<Q_t>(eQR.matrixQ() * Eigen::MatrixXd::Identity(m, n));

    auto & R1 = eQR.matrixR().template triangularView<Eigen::Upper>();
    Rmat_ = std::make_shared<R_type>(R1.block(0,0,n,n));
  }

  template <typename vector_t>
  void solveImpl(const vector_t & rhs, vector_t & y){
    ::rompp::qr::impl::solve(rhs, *Rmat_, y);
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
