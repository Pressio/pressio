
#ifndef QR_HOUSEHOLDER_EIGEN_DENSE_MATRIX_OUT_OF_PLACE_IMPL_HPP_
#define QR_HOUSEHOLDER_EIGEN_DENSE_MATRIX_OUT_OF_PLACE_IMPL_HPP_

#include <Eigen/QR>
#include "../../../CORE_OPS"
#include "../qr_rfactor_solve_impl.hpp"

namespace rompp{ namespace qr{ namespace impl{

/* partially specialize for when n and m are dynamic.
 * This does not mean that the wrapped matrix is dynamic, this just means
 * that m and n are passed as dynamic, so the Q factor will be a dynamic matrix.
*/
template< typename matrix_type, typename R_t, template <typename...> class Q_type>
class QRHouseholderDenseEigenMatrixWrapper<
  matrix_type, R_t, core::constants::dynamic, core::constants::dynamic, Q_type,
  core::meta::enable_if_t<
    core::meta::is_dense_matrix_wrapper_eigen<matrix_type>::value
    >
  >{

  using sc_t	     = typename core::details::traits<matrix_type>::scalar_t;
  using nat_mat_t    = typename core::details::traits<matrix_type>::wrapped_t;
  using factorizer_t = Eigen::HouseholderQR<nat_mat_t>;
  using Q_nat_t	     = Eigen::Matrix<sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  using Q_t	     = Q_type<Q_nat_t>;

  mutable std::shared_ptr<Q_t> Qmat_	     = {};
  mutable std::shared_ptr<factorizer_t> fct_ = {};

public:
  QRHouseholderDenseEigenMatrixWrapper() = default;
  ~QRHouseholderDenseEigenMatrixWrapper() = default;

  void computeThinOutOfPlace(matrix_type & A){
    auto rows = A.rows();
    auto cols = A.cols();
    fct_ = std::make_shared<factorizer_t>(*A.data());

    if (!Qmat_ or (Qmat_->data()->rows()!=rows and Qmat_->data()->cols()!=cols ) )
      Qmat_ = std::make_shared<Q_t>(rows,cols);

    *Qmat_->data() = fct_->householderQ() * Q_nat_t::Identity(rows,cols);
  }

  template < typename vector_in_t, typename vector_out_t>
  void project(const vector_in_t & vecIn,
  	       vector_out_t & vecOut) const{
    core::ops::dot( *this->Qmat_, vecIn, vecOut );
  }

  // non-type template n not used here
  template <typename vector_t>
  void doLinSolve(const vector_t & rhs, vector_t & y)const{
    auto vecSize = y.size();
    auto & Rm = fct_->matrixQR().block(0,0,vecSize,vecSize).
      template triangularView<Eigen::Upper>();
    *y.data() = Rm.solve(*rhs.data());
  }

  const Q_t & getCRefQFactor() const {
    return *this->Qmat_;
  }
};




/* partially specialize for when n and m are known at compile time.
 * If m,n known, just tells me that Q can be statically allocated.
 * It does NOT necessarily imply the wrapped matrix is static too,
 * the wrapped matrix can be dynamic, it does not matter.
*/
template< typename matrix_type, typename R_t,
	  int n, int m, template <typename...> class Q_type>
class QRHouseholderDenseEigenMatrixWrapper<
  matrix_type, R_t, n, m, Q_type,
  core::meta::enable_if_t<
    core::meta::is_dense_matrix_wrapper_eigen<matrix_type>::value and n >=1 and m>=1>
    >{

  using sc_t	     = typename core::details::traits<matrix_type>::scalar_t;
  using nat_mat_t    = typename core::details::traits<matrix_type>::wrapped_t;
  using factorizer_t = Eigen::HouseholderQR<nat_mat_t>;
  using Q_nat_t	     = Eigen::Matrix<sc_t, m, n>;
  using Q_t	     = Q_type<Q_nat_t>;

  Q_nat_t eigIdentity_			     = Q_nat_t::Identity();
  mutable std::shared_ptr<Q_t> Qmat_	     = {};
  mutable std::shared_ptr<factorizer_t> fct_ = {};

public:
  QRHouseholderDenseEigenMatrixWrapper() = default;
  ~QRHouseholderDenseEigenMatrixWrapper() = default;

  void computeThinOutOfPlace(matrix_type & A){
    auto rows = A.rows();
    auto cols = A.cols();
    // if doing thin, then m,n need to satisfy these constraints
    assert(rows == m);
    assert(cols == n);

    fct_ = std::make_shared<factorizer_t>(*A.data());

    if (!Qmat_ or
	(Qmat_->data()->rows()!=rows and
	 Qmat_->data()->cols()!=cols)){
      Qmat_ = std::make_shared<Q_t>(); // static allocation
    }

    *Qmat_->data() = fct_->householderQ() * eigIdentity_;
  }

  template < typename vector_in_t, typename vector_out_t>
  void project(const vector_in_t & vecIn,
  	       vector_out_t & vecOut) const{
    core::ops::dot( *this->Qmat_, vecIn, vecOut );
  }

  // non-type template n not used here
  template <typename vector_t>
  void doLinSolve(const vector_t & rhs, vector_t & y)const{
    auto vecSize = y.size();
    auto & Rm = fct_->matrixQR().block(0,0,vecSize,vecSize).
      template triangularView<Eigen::Upper>();
    *y.data() = Rm.solve(*rhs.data());
  }

  const Q_t & getCRefQFactor() const {
    return *this->Qmat_;
  }
};



}}} // end namespace rompp::qr::impl
#endif
