
#ifndef QR_IN_PLACE_HPP_
#define QR_IN_PLACE_HPP_

#include "./base/qr_in_place_base.hpp"
#include "./base/qr_solve_base.hpp"
#include "./base/qr_r_factor_base.hpp"
#include "qr_traits.hpp"
#include "qr_rfactor_solve_impl.hpp"

namespace pressio{ namespace qr{ namespace impl{

/* overload for R_type == void */
template<
  typename matrix_type, typename algo, int m,
  int n, template <typename...> class Q_type
  >
class QRSolver<
  matrix_type, algo, true, m, n, void, Q_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper<matrix_type>::value or
    containers::meta::is_matrix_wrapper<matrix_type>::value
    >
  > : public details::traits< QRSolver<matrix_type, algo,
				       true, m, n, void, Q_type>>::base_compute_t,
      public details::traits< QRSolver<matrix_type, algo,
				       true, m, n, void, Q_type>>::base_solve_t
{

  using this_t		= QRSolver<matrix_type, algo, true, m, n, void, Q_type>;
  using traits_t	= details::traits<this_t>;
  using base_compute_t	= typename traits_t::base_compute_t;
  using base_solve_t	= typename traits_t::base_solve_t;

  using impl_t	       = typename traits_t::impl_t;
  impl_t myImpl_;

  void computeThinImpl(matrix_type & A){
    myImpl_.computeThinInPlace(A);
  }

  template <typename vector_t>
  void solveImpl(const vector_t & rhs, vector_t & y) const {
    myImpl_.template doLinSolve<vector_t>(rhs, y);
  }

public:
  QRSolver() = default;
  ~QRSolver() = default;

private:
  friend base_compute_t;
  friend base_solve_t;
};



// /* overload for R_type != void */
// template<
//   typename matrix_type, typename algo, int m,
//   int n, typename R_type, template <typename...> class Q_type
//   >
// class QRSolver<
//   matrix_type, algo, true, m, n, R_type, Q_type,
//   ::pressio::mpl::enable_if_t<
//     meta::is_legitimate_r_type<R_type>::value and
//     containers::meta::is_multi_vector_wrapper<matrix_type>::value
//     >
//   > : public details::traits< QRSolver<matrix_type, algo,
// 				       true, m, n, R_type, Q_type>>::base_compute_t,
//       public details::traits< QRSolver<matrix_type, algo,
// 				       true, m, n, R_type, Q_type>>::base_solve_t,
//       public details::traits< QRSolver<matrix_type, algo,
// 				       true, m, n, R_type, Q_type>>::base_Rfactor_t
// {

//   using this_t	       = QRSolver<matrix_type, algo, true, m, n, R_type, Q_type>;
//   using traits_t       = details::traits<this_t>;
//   using base_compute_t = typename traits_t::base_compute_t;
//   using base_solve_t   = typename traits_t::base_solve_t;
//   using base_Rfactor_t = typename traits_t::base_Rfactor_t;

//   using impl_t	       = typename traits_t::impl_t;
//   impl_t myImpl_       = {};

//   void computeThinImpl(matrix_type & A){
//     myImpl_.computeThinInPlace(A);
//   }

//   // here we have the Rfactor, maybe we should call the solver
//   // native to the Rfactor??? think
//   template <typename vector_t>
//   void solveImpl(const vector_t & rhs, vector_t & y) const{
//     myImpl_.template doLinSolve<vector_t, n>(rhs, y);
//   }

//   const R_type & cRefRFactorImpl() const {
//     return myImpl_.getCRefRFactor();
//   }

// public:
//   QRSolver() = default;
//   ~QRSolver() = default;

// private:
//   friend base_compute_t;
//   friend base_Rfactor_t;
//   friend base_solve_t;
// };


}}} // end namespace pressio::qr::impl
#endif
