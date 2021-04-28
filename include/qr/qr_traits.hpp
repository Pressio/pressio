/*
//@HEADER
// ************************************************************************
//
// qr_traits.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef QR_QR_TRAITS_HPP_
#define QR_QR_TRAITS_HPP_

namespace pressio{ namespace qr{ namespace details{

/* common to all cases */
template<typename matrix_type, typename algo, bool in_place>
struct traits_shared_all{
  using matrix_t  = matrix_type;
  using algo_t	  = algo;
  using nat_mat_t = typename containers::details::traits<matrix_type>::wrapped_t;
  using sc_t	  = typename containers::details::traits<matrix_type>::scalar_t;
  static constexpr bool in_place_ = in_place;
};


/* helpers to define the type of the implementation class to use */
template <
  typename matrix_t,
  typename algo_tag,
  typename R_t,
  typename wrap_Q_type,
  template <typename...> class Q_type,
  typename enable = void>
struct impl_class_helper{};


#if defined PRESSIO_ENABLE_TPL_TRILINOS
template <
  typename matrix_t, typename R_t, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<
  matrix_t, qr::TSQR, R_t, wrap_Q_type, Q_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_multi_vector_wrapper_epetra<matrix_t>::value
    >
  >
{
  using impl_t = impl::EpetraMVTSQR<matrix_t, R_t, wrap_Q_type, Q_type>;
};

template <
  typename matrix_t, typename R_t, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<
  matrix_t, qr::TSQR, R_t, wrap_Q_type, Q_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_multi_vector_wrapper_tpetra<matrix_t>::value
    >
  >
{
  using impl_t = impl::TpetraMVTSQR<matrix_t, R_t, wrap_Q_type, Q_type>;
};

template <
  typename matrix_t, typename R_t, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<
  matrix_t, qr::TSQR, R_t, wrap_Q_type, Q_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_multi_vector_wrapper_tpetra_block<matrix_t>::value
    >
  >
{
  using impl_t = impl::TpetraBlockMVTSQR<matrix_t, R_t, wrap_Q_type, Q_type>;
};

template <
  typename matrix_t, typename R_t, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<
  matrix_t, qr::ModifiedGramSchmidt, R_t, wrap_Q_type, Q_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_multi_vector_wrapper_epetra<matrix_t>::value
    >
  >
{
  using impl_t = impl::ModGramSchmidtMVEpetra<matrix_t, R_t, wrap_Q_type, Q_type>;
};

template <
  typename matrix_t, typename R_t, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<
  matrix_t, qr::ModifiedGramSchmidt, R_t, wrap_Q_type, Q_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_multi_vector_wrapper_tpetra<matrix_t>::value
    >
  >
{
  using impl_t = impl::ModGramSchmidtMVTpetra<matrix_t, R_t, wrap_Q_type, Q_type>;
};


template <
  typename matrix_t, typename R_t, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<
  matrix_t, qr::Householder, R_t, wrap_Q_type, Q_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_multi_vector_wrapper_epetra<matrix_t>::value
    >
  >
{
  using impl_t = impl::EpetraMVHouseholderUsingEigen<matrix_t, R_t, Q_type>;
};


template <
  typename matrix_t, typename R_t, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<
  matrix_t, qr::Householder, R_t, wrap_Q_type, Q_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_multi_vector_wrapper_tpetra<matrix_t>::value
    >
  >{
  using impl_t = impl::TpetraMVHouseholderUsingEigen<matrix_t, R_t, Q_type>;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <
  typename matrix_t, typename R_t, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<
  matrix_t, qr::Householder, R_t, wrap_Q_type, Q_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_dense_matrix_wrapper_eigen<matrix_t>::value
    >
  >
{
  using impl_t = impl::QRHouseholderDenseEigenMatrixWrapper<matrix_t, R_t, Q_type>;
};

template <
  typename matrix_t, typename R_t, typename wrap_Q_type, template <typename...> class Q_type>
struct impl_class_helper<
  matrix_t, qr::Householder, R_t, wrap_Q_type, Q_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_multi_vector_wrapper_eigen<matrix_t>::value
    >
  >
{
  using impl_t = impl::QRHouseholderEigenMultiVectorWrapper<matrix_t, R_t, Q_type>;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_EIGEN
/*
 * specialize for:
 *	Eigen::DenseMatrixWrapper,
 *	R_type = void
 */
template<
  typename matrix_type,
  typename algo_t,
  bool in_place,
  template <typename...> class Q_type
  >
struct traits<
  impl::QRSolver<
    matrix_type, algo_t, in_place, void, Q_type>,
    ::pressio::mpl::enable_if_t<
      containers::predicates::is_dense_matrix_wrapper_eigen<matrix_type>::value or
      containers::predicates::is_multi_vector_wrapper_eigen<matrix_type>::value
      >
  > : traits_shared_all<matrix_type, algo_t, in_place>
{

  static_assert( std::is_same<algo_t, ModifiedGramSchmidt>::value or
		 std::is_same<algo_t, Householder>::value,
  "Currently, only ModifiedGramSchmidt and Householder are available for Eigen matrices.");

  using traits_all_t	= traits_shared_all<matrix_type, algo_t, in_place>;
  using typename traits_all_t::matrix_t;
  using typename traits_all_t::sc_t;
  using typename traits_all_t::nat_mat_t;

  using nat_Q_t		= Eigen::Matrix<sc_t, Eigen::Dynamic, Eigen::Dynamic>;
  using Q_t		= Q_type<nat_Q_t>;

  using concrete_t	= impl::QRSolver<matrix_type, algo_t, in_place, void, Q_type>;
  using inplace_base_t  = QRInPlaceBase<concrete_t, matrix_type>;
  using outplace_base_t = QROutOfPlaceBase<concrete_t, matrix_type, Q_t>;
  using base_compute_t	= typename std::conditional<in_place, inplace_base_t, outplace_base_t>::type;
  using base_solve_t	= QRSolveBase<concrete_t>;
  using impl_t		= typename impl_class_helper<matrix_t, algo_t, void, nat_Q_t, Q_type>::impl_t;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS

/*
 * traits_shared_trilinos_mv
 */
template<
  typename matrix_type, template <typename...> class Q_type,
  ::pressio::mpl::enable_if_t<
    containers::predicates::is_multi_vector_wrapper_epetra<matrix_type>::value or
    containers::predicates::is_multi_vector_wrapper_tpetra<matrix_type>::value or
    containers::predicates::is_multi_vector_wrapper_tpetra_block<matrix_type>::value, int > = 0
  >
struct traits_shared_trilinos_mv{
  using MV_t	 = typename containers::details::traits<matrix_type>::wrapped_t;
  using LO_t	 = typename containers::details::traits<matrix_type>::local_ordinal_t;
  using GO_t	 = typename containers::details::traits<matrix_type>::global_ordinal_t;
  using map_t	 = typename containers::details::traits<matrix_type>::data_map_t;
  using Q_t	 = Q_type<MV_t>;
};


/*
 * specialize for Epetra::MultiVector, R_type = void
 */
template<typename matrix_type, typename algo_t, bool in_place, template <typename...> class Q_type>
struct traits<
  impl::QRSolver<
    matrix_type, algo_t, in_place, void, Q_type>,
    ::pressio::mpl::enable_if_t<
      containers::predicates::is_multi_vector_wrapper_epetra<matrix_type>::value
      >
  > : traits_shared_all<matrix_type, algo_t, in_place>,
  traits_shared_trilinos_mv<matrix_type, Q_type>
{

  static_assert( std::is_same<algo_t, ModifiedGramSchmidt>::value or
		 std::is_same<algo_t, Householder>::value or
		 std::is_same<algo_t, TSQR>::value,
		 "Currently, only TSQR, ModifiedGramSchmidt and Householder are available for Epetra dense matrices. Use TSQR because it is fast and accurate. ModifiedGramSchmidt and Householder are just here for testing purposes. ");

  using traits_all_t  = traits_shared_all<matrix_type, algo_t, in_place>;
  using traits_tril_t = traits_shared_trilinos_mv<matrix_type, Q_type>;
  using typename traits_all_t::matrix_t;
  using typename traits_all_t::sc_t;
  using typename traits_tril_t::Q_t;
  using typename traits_tril_t::MV_t;

  using concrete_t	= impl::QRSolver<matrix_type, algo_t, in_place, void, Q_type>;
  using inplace_base_t  = QRInPlaceBase<concrete_t, matrix_type>;
  using outplace_base_t = QROutOfPlaceBase<concrete_t, matrix_type, Q_t>;
  using base_compute_t	= typename std::conditional<in_place, inplace_base_t, outplace_base_t>::type;
  using base_solve_t	= QRSolveBase<concrete_t>;
  using impl_t		= typename impl_class_helper<matrix_t, algo_t, void, MV_t, Q_type>::impl_t;
};


/*
 * specialize for Tpetra::MultiVector, R_type = void
 */
template<typename matrix_type, typename algo_t, bool in_place, template <typename...> class Q_type>
struct traits<
  impl::QRSolver<
    matrix_type, algo_t, in_place, void, Q_type>,
    ::pressio::mpl::enable_if_t<
      containers::predicates::is_multi_vector_wrapper_tpetra<matrix_type>::value
      >
  > : traits_shared_all<matrix_type, algo_t, in_place>,
  traits_shared_trilinos_mv<matrix_type, Q_type>
{

  static_assert( std::is_same<algo_t, ModifiedGramSchmidt>::value or
		 std::is_same<algo_t, Householder>::value or
		 std::is_same<algo_t, TSQR>::value,
		 "Currently, only TSQR, ModifiedGramSchmidt and Householder are available for Tpetra dense matrices. Use TSQR because it is fast and accurate. ModifiedGramSchmidt and Householder are just here for testing purposes. ");

  using traits_all_t  = traits_shared_all<matrix_type, algo_t, in_place>;
  using traits_tril_t = traits_shared_trilinos_mv<matrix_type, Q_type>;

  using typename traits_all_t::matrix_t;
  using typename traits_all_t::sc_t;
  using typename traits_tril_t::Q_t;
  using typename traits_tril_t::MV_t;
  using node_t = typename containers::details::traits<matrix_type>::node_t;
  using hexsp  = typename containers::details::traits<matrix_type>::host_exec_space_t;

  using concrete_t	= impl::QRSolver<matrix_type, algo_t, in_place, void, Q_type>;
  using inplace_base_t  = QRInPlaceBase<concrete_t, matrix_type>;
  using outplace_base_t = QROutOfPlaceBase<concrete_t, matrix_type, Q_t>;

  using base_compute_t	= typename std::conditional<in_place, inplace_base_t, outplace_base_t>::type;
  using base_solve_t	= QRSolveBase<concrete_t>;
  using impl_t		= typename impl_class_helper<matrix_t, algo_t, void, MV_t, Q_type>::impl_t;
};


/*
 * specialize for Tpetra::BlockMultiVector, R_type = void
 */
template<typename matrix_type, typename algo_t, bool in_place, template <typename...> class Q_type>
struct traits<
  impl::QRSolver<
    matrix_type, algo_t, in_place, void, Q_type>,
    ::pressio::mpl::enable_if_t<
      containers::predicates::is_multi_vector_wrapper_tpetra_block<matrix_type>::value
      >
  > : traits_shared_all<matrix_type, algo_t, in_place>,
  traits_shared_trilinos_mv<matrix_type, Q_type>
{

  static_assert( std::is_same<algo_t, TSQR>::value,
		 "Currently, only TSQR is available for BlockTpetra dense matrices.");

  using traits_all_t  = traits_shared_all<matrix_type, algo_t, in_place>;
  using traits_tril_t = traits_shared_trilinos_mv<matrix_type, Q_type>;

  using typename traits_all_t::matrix_t;
  using typename traits_all_t::sc_t;
  using typename traits_tril_t::Q_t;
  using typename traits_tril_t::MV_t;
  using node_t = typename containers::details::traits<matrix_type>::node_t;
  using hexsp  = typename containers::details::traits<matrix_type>::host_exec_space_t;

  using concrete_t	= impl::QRSolver<matrix_type, algo_t, in_place, void, Q_type>;
  using inplace_base_t  = QRInPlaceBase<concrete_t, matrix_type>;
  using outplace_base_t = QROutOfPlaceBase<concrete_t, matrix_type, Q_t>;

  using base_compute_t	= typename std::conditional<in_place,inplace_base_t, outplace_base_t>::type;
  using base_solve_t	= QRSolveBase<concrete_t>;
  using impl_t		= typename impl_class_helper<matrix_t, algo_t, void, MV_t, Q_type>::impl_t;
};
#endif //PRESSIO_ENABLE_TPL_TRILINOS

}}}//end namespace pressio::qr::details






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
//     ::pressio::mpl::enable_if_t<
//       containers::predicates::is_multi_vector_wrapper_epetra<matrix_type>::value and
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
//     ::pressio::mpl::enable_if_t<
//       containers::predicates::is_multi_vector_wrapper_tpetra<matrix_type>::value and
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
#endif  // QR_QR_TRAITS_HPP_
