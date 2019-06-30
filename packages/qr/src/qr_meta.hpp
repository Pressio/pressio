
#ifndef QR_META_HPP_
#define QR_META_HPP_

#include "qr_ConfigDefs.hpp"
#include "qr_algorithms_tags.hpp"
#include "../../containers/src/multi_vector/containers_multi_vector_meta.hpp"
#include "../../containers/src/vector/containers_vector_meta.hpp"
#include "../../containers/src/matrix/containers_matrix_meta.hpp"


namespace rompp{ namespace qr{ namespace meta {

template <typename T, typename enable = void>
struct is_legitimate_r_type : std::false_type {};

template <typename T>
struct is_legitimate_r_type<T,
	 ::rompp::mpl::enable_if_t<
	   containers::meta::is_matrix_wrapper<T>::value and
	   containers::details::traits<T>::is_shared_mem and
	   containers::details::traits<T>::is_dense
	   >
      > : std::true_type{};



template <typename T, typename Q_T, typename enable = void>
struct is_legitimate_vector_type_for_qr_project : std::false_type {};

template <typename T, typename Q_t>
struct is_legitimate_vector_type_for_qr_project<T, Q_t,
	 ::rompp::mpl::enable_if_t<
	   containers::meta::is_vector_wrapper<T>::value and
	   // the vector type should be from same package as Q
	   containers::details::traits<T>::wrapped_package_identifier ==
	   containers::details::traits<Q_t>::wrapped_package_identifier
	 >
      > : std::true_type{};


#if defined HAVE_TRILINOS
template <typename algo_t, typename enable = void>
struct is_legitimate_algo_for_epetra_mv : std::false_type {};

template <typename algo_t>
struct is_legitimate_algo_for_epetra_mv<algo_t,
	 ::rompp::mpl::enable_if_t<
	   std::is_same<algo_t, ::rompp::qr::Householder>::value
	   or std::is_same<algo_t, ::rompp::qr::TSQR>::value
	 >
      > : std::true_type{};
#endif


#if defined HAVE_TRILINOS
template <typename algo_t, typename enable = void>
struct is_legitimate_algo_for_tpetra_mv : std::false_type {};

template <typename algo_t>
struct is_legitimate_algo_for_tpetra_mv<algo_t,
	 ::rompp::mpl::enable_if_t<
	   std::is_same<algo_t, ::rompp::qr::Householder>::value
	   or std::is_same<algo_t, ::rompp::qr::TSQR>::value
	   >
        > : std::true_type{};
#endif

}}}//end namespace rompp::qr::meta
#endif
