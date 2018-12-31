
#ifndef QR_EIGEN_DENSE_MATRIX_OUT_OF_PLACE_HPP_
#define QR_EIGEN_DENSE_MATRIX_OUT_OF_PLACE_HPP_

#include "./base/qr_out_of_place_base.hpp"
#include "./base/qr_solve_base.hpp"
#include "./base/qr_r_factor_base.hpp"
#include "qr_traits.hpp"
#include "qr_rfactor_solve_impl.hpp"
#include <Eigen/QR>
#include "../../CORE_OPS"

namespace rompp{ namespace qr{ namespace impl{

/*
 *   overload for R_type == void, in_place = false
 */
template<
  typename matrix_type, typename algo, int m,
  int n, template <typename...> class Q_type
  >
class QRSolver<
  matrix_type, algo, false, m, n, void, Q_type,
  core::meta::enable_if_t<
    core::meta::is_eigen_dense_matrix_wrapper<matrix_type>::value and
    std::is_same<algo, ::rompp::qr::Householder>::value
    >
  > : public details::traits< QRSolver<matrix_type, algo,
				       false, m, n, void, Q_type>>::base_compute_t,
      public details::traits< QRSolver<matrix_type, algo,
				       false, m, n, void, Q_type>>::base_solve_t
{

  using this_t	       = QRSolver<matrix_type, algo, false, m, n, void, Q_type>;
  using traits_t       = details::traits<this_t>;
  using base_compute_t = typename traits_t::base_compute_t;
  using base_solve_t   = typename traits_t::base_solve_t;
  using Q_t	       = typename traits_t::Q_t;
  using nat_mat_t      = typename traits_t::nat_mat_t;
  using factorizer_t   = Eigen::HouseholderQR<nat_mat_t>;

  mutable std::shared_ptr<Q_t> Qmat_	     = {};
  mutable std::shared_ptr<factorizer_t> fct_ = {};

public:
  QRSolver() = default;
  ~QRSolver() = default;

private:
  void computeThinImpl(matrix_type & A){
    auto rows = A.rows();
    auto cols = A.cols();

    fct_ = std::make_shared<factorizer_t>(*A.data());

    if (!Qmat_ or (Qmat_->data()->rows()!=rows and Qmat_->data()->cols()!=cols ) )
      Qmat_ = std::make_shared<Q_t>(rows,cols);

    *Qmat_->data() = fct_->householderQ() * Eigen::MatrixXd::Identity(rows,cols);
  }

  template < typename vector_in_t, typename vector_out_t>
  void projectImpl(const vector_in_t & vecIn,
  		   vector_out_t & vecOut) const{
    core::ops::dot( *this->Qmat_, vecIn, vecOut );
  }

  template <typename vector_t>
  void solveImpl(const vector_t & rhs, vector_t & y)const{
    auto vecSize = y.size();
    auto & Rm = fct_->matrixQR().block(0,0,vecSize,vecSize).template triangularView<Eigen::Upper>();
    *y.data() = Rm.solve(*rhs.data());
  }

  const Q_t & cRefQFactorImpl() const {
    return *this->Qmat_;
  }

private:
  friend base_compute_t;
  friend base_solve_t;

};


}}} // end namespace rompp::qr::impl
#endif
