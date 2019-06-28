
#ifndef QR_TRAITS_HPP_
#define QR_TRAITS_HPP_

#include "qr_forward_declarations.hpp"
#include "qr_meta.hpp"

namespace rompp{ namespace qr{ namespace details{

/* common to all cases */
template< typename matrix_type, typename algo,
	  bool in_place, int m, int n >
struct traits_shared_all{

  using matrix_t  = matrix_type;
  using algo_t	  = algo;
  using nat_mat_t = typename algebra::details::traits<matrix_type>::wrapped_t;
  using sc_t	  = typename algebra::details::traits<matrix_type>::scalar_t;

  static constexpr bool in_place_ = in_place;
  static constexpr int m_	  = m;
  static constexpr int n_	  = n;
};


/* helpers to define the type of the implementation class to use */
template <typename matrix_t, typename algo_tag, typename R_t,
	  int n, int m, typename wrap_Q_type, template <typename...> class Q_type,
	  typename enable = void>
struct impl_class_helper{};


#if defined HAVE_TRILINOS
template <typename matrix_t, typename R_t,
	  int n, int m, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<matrix_t, qr::TSQR, R_t, n, m, wrap_Q_type, Q_type,
			 ::rompp::mpl::enable_if_t<
			   algebra::meta::is_multi_vector_wrapper_epetra<matrix_t>::value
			   >>{
  using impl_t = impl::EpetraMVTSQR<matrix_t, R_t, n, m, wrap_Q_type, Q_type>;
};

template <typename matrix_t, typename R_t,
	  int n, int m, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<matrix_t, qr::TSQR, R_t, n, m, wrap_Q_type, Q_type,
			 ::rompp::mpl::enable_if_t<
			   algebra::meta::is_multi_vector_wrapper_tpetra<matrix_t>::value
			   >>{
  using impl_t = impl::TpetraMVTSQR<matrix_t, R_t, n, m, wrap_Q_type, Q_type>;
};

template <typename matrix_t, typename R_t,
	  int n, int m, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<matrix_t, qr::TSQR, R_t, n, m, wrap_Q_type, Q_type,
			 ::rompp::mpl::enable_if_t<
			   algebra::meta::is_multi_vector_wrapper_tpetra_block<matrix_t>::value
			   >>{
  using impl_t = impl::TpetraBlockMVTSQR<matrix_t, R_t, n, m, wrap_Q_type, Q_type>;
};

template <typename matrix_t, typename R_t,
	  int n, int m, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<matrix_t, qr::ModifiedGramSchmidt, R_t, n, m, wrap_Q_type, Q_type,
			 ::rompp::mpl::enable_if_t<
			   algebra::meta::is_multi_vector_wrapper_epetra<matrix_t>::value
			   >>{
  using impl_t = impl::ModGramSchmidtMVEpetra<matrix_t, R_t, n, m, wrap_Q_type, Q_type>;
};

template <typename matrix_t, typename R_t,
	  int n, int m, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<matrix_t, qr::ModifiedGramSchmidt, R_t, n, m, wrap_Q_type, Q_type,
			 ::rompp::mpl::enable_if_t<
			   algebra::meta::is_multi_vector_wrapper_tpetra<matrix_t>::value
			   >>{
  using impl_t = impl::ModGramSchmidtMVTpetra<matrix_t, R_t, n, m, wrap_Q_type, Q_type>;
};


template <typename matrix_t, typename R_t,
	  int n, int m, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<matrix_t, qr::Householder, R_t, n, m, wrap_Q_type, Q_type,
			 ::rompp::mpl::enable_if_t<
			   algebra::meta::is_multi_vector_wrapper_epetra<matrix_t>::value
			   >>{
  using impl_t = impl::EpetraMVHouseholderUsingEigen<matrix_t, R_t, n, m, Q_type>;
};


template <typename matrix_t, typename R_t,
	  int n, int m, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<matrix_t, qr::Householder, R_t, n, m, wrap_Q_type, Q_type,
			 ::rompp::mpl::enable_if_t<
			   algebra::meta::is_multi_vector_wrapper_tpetra<matrix_t>::value
			   >>{
  using impl_t = impl::TpetraMVHouseholderUsingEigen<matrix_t, R_t, n, m, Q_type>;
};
#endif


template <typename matrix_t, typename R_t,
	  int n, int m, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<matrix_t, qr::Householder, R_t, n, m, wrap_Q_type, Q_type,
			 ::rompp::mpl::enable_if_t<
			   algebra::meta::is_dense_matrix_wrapper_eigen<matrix_t>::value
			   >>{
  using impl_t = impl::QRHouseholderDenseEigenMatrixWrapper<matrix_t, R_t, n, m, Q_type>;
};


/*
 * specialize for:
 *	Eigen::DenseMatrixWrapper,
 *	R_type = void
 */
template<
  typename matrix_type, typename algo_t, bool in_place, int m, int n,
  template <typename...> class Q_type
  >
struct traits<
  impl::QRSolver<
    matrix_type, algo_t, in_place, m, n, void, Q_type>,
    ::rompp::mpl::enable_if_t<
      algebra::meta::is_dense_matrix_wrapper_eigen<matrix_type>::value
      >
  > : traits_shared_all<matrix_type, algo_t, in_place, m, n>{

  static_assert( std::is_same<algo_t, ModifiedGramSchmidt>::value or
		 std::is_same<algo_t, Householder>::value,
  "Currently, only ModifiedGramSchmidt and Householder are available for Eigen matrices.");

  using traits_all_t	= traits_shared_all<matrix_type, algo_t, in_place, m, n>;
  using typename traits_all_t::matrix_t;
  using typename traits_all_t::sc_t;
  using typename traits_all_t::nat_mat_t;

  using nat_Q_t		= Eigen::Matrix<sc_t, m, n>;
  using Q_t		= Q_type<nat_Q_t>;

  using concrete_t	= impl::QRSolver<matrix_type, algo_t,
					in_place, m, n, void, Q_type>;
  using inplace_base_t  = QRInPlaceBase<concrete_t, matrix_type>;
  using outplace_base_t = QROutOfPlaceBase<concrete_t, matrix_type, Q_t>;
  using base_compute_t	= typename std::conditional<in_place,
						    inplace_base_t,
						    outplace_base_t>::type;
  using base_solve_t	= QRSolveBase<concrete_t>;
  using impl_t		= typename impl_class_helper<matrix_t, algo_t, void, n, m,
						     nat_Q_t, Q_type>::impl_t;
};


#ifdef HAVE_TRILINOS

/*
 * traits_shared_trilinos_mv
 */
template<
  typename matrix_type, template <typename...> class Q_type,
  ::rompp::mpl::enable_if_t<
    algebra::meta::is_multi_vector_wrapper_epetra<matrix_type>::value or
    algebra::meta::is_multi_vector_wrapper_tpetra<matrix_type>::value or
    algebra::meta::is_multi_vector_wrapper_tpetra_block<matrix_type>::value
    > * = nullptr
  >
struct traits_shared_trilinos_mv{

  using MV_t	 = typename algebra::details::traits<matrix_type>::wrapped_t;
  using LO_t	 = typename algebra::details::traits<matrix_type>::local_ordinal_t;
  using GO_t	 = typename algebra::details::traits<matrix_type>::global_ordinal_t;
  using map_t	 = typename algebra::details::traits<matrix_type>::data_map_t;
  using Q_t	 = Q_type<MV_t>;
};


/*
 * specialize for Epetra::MultiVector, R_type = void
 */
template<
  typename matrix_type, typename algo_t, bool in_place, int m,
  int n, template <typename...> class Q_type
  >
struct traits<
  impl::QRSolver<
    matrix_type, algo_t, in_place, m, n, void, Q_type>,
    ::rompp::mpl::enable_if_t<
      algebra::meta::is_multi_vector_wrapper_epetra<matrix_type>::value
      >
  > : traits_shared_all<matrix_type, algo_t, in_place, m, n>,
  traits_shared_trilinos_mv<matrix_type, Q_type>{

  static_assert( std::is_same<algo_t, ModifiedGramSchmidt>::value or
		 std::is_same<algo_t, Householder>::value or
		 std::is_same<algo_t, TSQR>::value,
		 "Currently, only TSQR, ModifiedGramSchmidt and Householder are available for Epetra dense matrices. Use TSQR because it is fast and accurate. ModifiedGramSchmidt and Householder are just here for testing purposes. ");

  using traits_all_t  = traits_shared_all<matrix_type, algo_t, in_place, m, n>;
  using traits_tril_t = traits_shared_trilinos_mv<matrix_type, Q_type>;
  using typename traits_all_t::matrix_t;
  using typename traits_all_t::sc_t;
  using typename traits_tril_t::Q_t;
  using typename traits_tril_t::MV_t;

  using concrete_t	= impl::QRSolver<matrix_type, algo_t,
					in_place, m, n, void, Q_type>;
  using inplace_base_t  = QRInPlaceBase<concrete_t, matrix_type>;
  using outplace_base_t = QROutOfPlaceBase<concrete_t, matrix_type, Q_t>;
  using base_compute_t	= typename std::conditional<in_place,
						    inplace_base_t,
						    outplace_base_t>::type;
  using base_solve_t	= QRSolveBase<concrete_t>;
  using impl_t		= typename impl_class_helper<matrix_t, algo_t, void, n, m,
						     MV_t, Q_type>::impl_t;
};


/*
 * specialize for Tpetra::MultiVector, R_type = void
 */
template<
  typename matrix_type, typename algo_t, bool in_place, int m,
  int n, template <typename...> class Q_type
  >
struct traits<
  impl::QRSolver<
    matrix_type, algo_t, in_place, m, n, void, Q_type>,
    ::rompp::mpl::enable_if_t<
      algebra::meta::is_multi_vector_wrapper_tpetra<matrix_type>::value
      >
  > : traits_shared_all<matrix_type, algo_t, in_place, m, n>,
  traits_shared_trilinos_mv<matrix_type, Q_type>{


  static_assert( std::is_same<algo_t, ModifiedGramSchmidt>::value or
		 std::is_same<algo_t, Householder>::value or
		 std::is_same<algo_t, TSQR>::value,
		 "Currently, only TSQR, ModifiedGramSchmidt and Householder are available for Tpetra dense matrices. Use TSQR because it is fast and accurate. ModifiedGramSchmidt and Householder are just here for testing purposes. ");

  using traits_all_t  = traits_shared_all<matrix_type, algo_t, in_place, m, n>;
  using traits_tril_t = traits_shared_trilinos_mv<matrix_type, Q_type>;

  using typename traits_all_t::matrix_t;
  using typename traits_all_t::sc_t;
  using typename traits_tril_t::Q_t;
  using typename traits_tril_t::MV_t;
  using node_t = typename algebra::details::traits<matrix_type>::node_t;
  using hexsp  = typename algebra::details::traits<matrix_type>::host_exec_space_t;

  using concrete_t	= impl::QRSolver<matrix_type, algo_t,
					in_place, m, n, void, Q_type>;
  using inplace_base_t  = QRInPlaceBase<concrete_t, matrix_type>;
  using outplace_base_t = QROutOfPlaceBase<concrete_t, matrix_type, Q_t>;

  using base_compute_t	= typename std::conditional<in_place,
						    inplace_base_t,
						    outplace_base_t>::type;
  using base_solve_t	= QRSolveBase<concrete_t>;
  using impl_t		= typename impl_class_helper<matrix_t, algo_t, void,
						     n, m, MV_t, Q_type>::impl_t;
};


/*
 * specialize for Tpetra::BlockMultiVector, R_type = void
 */
template<
  typename matrix_type,
  typename algo_t,
  bool in_place,
  int m, int n,
  template <typename...> class Q_type
  >
struct traits<
  impl::QRSolver<
    matrix_type, algo_t, in_place, m, n, void, Q_type>,
    ::rompp::mpl::enable_if_t<
      algebra::meta::is_multi_vector_wrapper_tpetra_block<matrix_type>::value
      >
  > : traits_shared_all<matrix_type, algo_t, in_place, m, n>,
  traits_shared_trilinos_mv<matrix_type, Q_type>{

  static_assert( std::is_same<algo_t, TSQR>::value,
		 "Currently, only TSQR is available for BlockTpetra dense matrices.");

  using traits_all_t  = traits_shared_all<matrix_type, algo_t, in_place, m, n>;
  using traits_tril_t = traits_shared_trilinos_mv<matrix_type, Q_type>;

  using typename traits_all_t::matrix_t;
  using typename traits_all_t::sc_t;
  using typename traits_tril_t::Q_t;
  using typename traits_tril_t::MV_t;
  using node_t = typename algebra::details::traits<matrix_type>::node_t;
  using hexsp  = typename algebra::details::traits<matrix_type>::host_exec_space_t;

  using concrete_t	= impl::QRSolver<matrix_type, algo_t,
					in_place, m, n, void, Q_type>;
  using inplace_base_t  = QRInPlaceBase<concrete_t, matrix_type>;
  using outplace_base_t = QROutOfPlaceBase<concrete_t, matrix_type, Q_t>;

  using base_compute_t	= typename std::conditional<in_place,
						    inplace_base_t,
						    outplace_base_t>::type;
  using base_solve_t	= QRSolveBase<concrete_t>;
  using impl_t		= typename impl_class_helper<matrix_t, algo_t, void,
						     n, m, MV_t, Q_type>::impl_t;
};


#endif //HAVE_TRILINOS

}}}//end namespace rompp::qr::details
#endif






// /*
//  * specialize for Epetra::MultiVector, R_type != void
//  */
// template<
//   typename matrix_type, typename algo_t, bool in_place, int m, int n,
//   typename R_type, template <typename...> class Q_type
//   >
// struct traits<
//   impl::QRSolver<
//     matrix_type, algo_t, in_place, m, n, R_type, Q_type>,
//     ::rompp::mpl::enable_if_t<
//       algebra::meta::is_multi_vector_wrapper_epetra<matrix_type>::value and
//       meta::is_legitimate_r_type<R_type>::value
//       >
//   > : traits_shared_all<matrix_type, algo_t, in_place, m, n>,
//   traits_shared_trilinos_mv<matrix_type, Q_type>{

//   using traits_all_t  = traits_shared_all<matrix_type, algo_t, in_place, m, n>;
//   using traits_tril_t = traits_shared_trilinos_mv<matrix_type, Q_type>;

//   using matrix_t	= typename traits_all_t::matrix_t;
//   using sc_t		= typename traits_all_t::sc_t;
//   using Q_t		= typename traits_tril_t::Q_t;
//   using MV_t		= typename traits_tril_t::MV_t;
//   using R_t		= R_type;

//   using concrete_t	= impl::QRSolver<matrix_type, algo_t,
// 					in_place, m, n, R_type, Q_type>;
//   using inplace_base_t  = QRInPlaceBase<concrete_t, matrix_type>;
//   using outplace_base_t = QROutOfPlaceBase<concrete_t, matrix_type, Q_t>;

//   using base_compute_t	= typename std::conditional<in_place,
// 						    inplace_base_t,
// 						    outplace_base_t>::type;
//   using base_solve_t	= QRSolveBase<concrete_t>;
//   using base_Rfactor_t	= RFactorBase<concrete_t, R_t>;
//   using impl_t		= typename impl_class_helper<matrix_t, algo_t, Q_t,
// 						     R_type, sc_t, MV_t>::impl_t;
// };



// /*
//  * specialize for Tpetra::MultiVector, R_type != void
//  */
// template<
//   typename matrix_type, typename algo_t, bool in_place, int m, int n,
//   typename R_type, template <typename...> class Q_type
//   >
// struct traits<
//   impl::QRSolver<
//     matrix_type, algo_t, in_place, m, n, R_type, Q_type>,
//     ::rompp::mpl::enable_if_t<
//       algebra::meta::is_multi_vector_wrapper_tpetra<matrix_type>::value and
//       meta::is_legitimate_r_type<R_type>::value
//       >
//   > : traits_shared_all<matrix_type, algo_t, in_place, m, n>,
//   traits_shared_trilinos_mv<matrix_type, Q_type>{

//   using traits_all_t  = traits_shared_all<matrix_type, algo_t, in_place, m, n>;
//   using traits_tril_t = traits_shared_trilinos_mv<matrix_type, Q_type>;

//   using matrix_t	= typename traits_all_t::matrix_t;
//   using sc_t		= typename traits_all_t::sc_t;
//   using Q_t		= typename traits_tril_t::Q_t;
//   using MV_t		= typename traits_tril_t::MV_t;
//   using R_t		= R_type;

//   using concrete_t	= impl::QRSolver<matrix_type, algo_t,
// 					in_place, m, n, R_type, Q_type>;
//   using inplace_base_t  = QRInPlaceBase<concrete_t, matrix_type>;
//   using outplace_base_t = QROutOfPlaceBase<concrete_t, matrix_type, Q_t>;

//   using base_compute_t	= typename std::conditional<in_place,
// 						    inplace_base_t,
// 						    outplace_base_t>::type;
//   using base_solve_t	= QRSolveBase<concrete_t>;
//   using base_Rfactor_t	= RFactorBase<concrete_t, R_t>;
//   using impl_t		= typename impl_class_helper<matrix_t, algo_t, Q_t,
// 						     R_type, sc_t, MV_t>::impl_t;
// };
