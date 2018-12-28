
#ifdef HAVE_TRILINOS
#ifndef QR_TRILINOS_MULTI_VECTOR_TSQR_OUT_OF_PLACE_HPP_
#define QR_TRILINOS_MULTI_VECTOR_TSQR_OUT_OF_PLACE_HPP_

#include "./base/qr_out_of_place_base.hpp"
#include "./base/qr_solve_base.hpp"
#include "./base/qr_r_factor_base.hpp"
#include "qr_traits.hpp"
#include "qr_rfactor_solve_impl.hpp"

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
    core::meta::is_epetra_multi_vector_wrapper<matrix_type>::value or
    core::meta::is_tpetra_multi_vector_wrapper<matrix_type>::value
    >
  > : public details::traits< QRSolver<matrix_type, algo,
				       false, m, n, void, Q_type>>::base_compute_t,
      public details::traits< QRSolver<matrix_type, algo,
				       false, m, n, void, Q_type>>::base_solve_t,
      public details::traits< QRSolver<matrix_type, algo,
				       false, m, n, void, Q_type>>::impl_t
{

  using this_t	       = QRSolver<matrix_type, algo, false, m, n, void, Q_type>;
  using traits_t       = details::traits<this_t>;
  using base_compute_t = typename traits_t::base_compute_t;
  using base_solve_t   = typename traits_t::base_solve_t;
  using Q_t	       = typename traits_t::Q_t;

public:
  QRSolver() = default;
  ~QRSolver() = default;

private:
  void computeThinImpl(matrix_type & A){
    this->computeThinOutOfPlace(A);
  }

  template < typename vector_in_t, typename vector_out_t>
  void projectImpl(const vector_in_t & vecIn,
  		   vector_out_t & vecOut) const{
    core::ops::dot( *this->Qmat_, vecIn, vecOut );
  }

  template <typename vector_t>
  void solveImpl(const vector_t & rhs, vector_t & y)const{
    this->template doLinSolve<vector_t, n>(rhs, y);
  }

  const Q_t & cRefQFactorImpl() const {
    return this->getCRefQFactor();
  }

private:
  friend base_compute_t;
  friend base_solve_t;

};



/*
 *   overload for R_type != void, in_place = false
 */
template<
  typename matrix_type, typename algo, int m,
  int n, typename R_type, template <typename...> class Q_type
  >
class QRSolver<
  matrix_type, algo, false, m, n, R_type, Q_type,
  core::meta::enable_if_t<
    meta::is_legitimate_r_type<R_type>::value and
    (core::meta::is_epetra_multi_vector_wrapper<matrix_type>::value or
     core::meta::is_tpetra_multi_vector_wrapper<matrix_type>::value)
    >
  > : public details::traits< QRSolver<matrix_type, algo,
				       false, m, n, R_type, Q_type>>::base_compute_t,
      public details::traits< QRSolver<matrix_type, algo,
				       false, m, n, R_type, Q_type>>::base_solve_t,
      public details::traits< QRSolver<matrix_type, algo,
				       false, m, n, R_type, Q_type>>::base_Rfactor_t,
      public details::traits< QRSolver<matrix_type, algo,
				       false, m, n, R_type, Q_type>>::impl_t
{

  using this_t	       = QRSolver<matrix_type, algo, false, m, n, R_type, Q_type>;
  using traits_t       = details::traits<this_t>;
  using base_compute_t = typename traits_t::base_compute_t;
  using base_solve_t   = typename traits_t::base_solve_t;
  using base_Rfactor_t = typename traits_t::base_Rfactor_t;
  using Q_t	       = typename traits_t::Q_t;

  void computeThinImpl(matrix_type & A){
    this->computeThinOutOfPlace(A);
  }

  const Q_t & cRefQFactorImpl() const {
    return this->getCRefQFactor();
  }

  template < typename vector_in_t, typename vector_out_t>
  void projectImpl(const vector_in_t & vecIn,
  		   vector_out_t & vecOut) const{
    core::ops::dot( *this->Qmat_, vecIn, vecOut );
  }

  // here we have the Rfactor, maybe we should call the solver
  // native to the Rfactor??? think
  template <typename vector_t>
  void solveImpl(const vector_t & rhs, vector_t & y)const {
    this->template doLinSolve<vector_t, n>(rhs, y);
  }

  const R_type & cRefRFactorImpl() const {
    return this->getCRefRFactor();
  }

public:
  QRSolver() = default;
  ~QRSolver() = default;

private:
  friend base_compute_t;
  friend base_Rfactor_t;
  friend base_solve_t;
};


}}} // end namespace rompp::qr::impl
#endif
#endif
