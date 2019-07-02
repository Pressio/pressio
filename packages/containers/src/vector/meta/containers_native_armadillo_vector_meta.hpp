
#ifdef HAVE_ARMADILLO
#ifndef CONTAINERS_NATIVE_ARMADILLO_VECTOR_META_HPP_
#define CONTAINERS_NATIVE_ARMADILLO_VECTOR_META_HPP_

#include "../meta/containers_meta_basic.hpp"
#include <armadillo>

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_column_vector_armadillo : std::false_type {};

template <typename T>
struct is_column_vector_armadillo<T,
	 ::pressio::mpl::enable_if_t<
	   std::is_same<T,
     	    arma::Col<typename T::elem_type>
			>::value
	   >
      > : std::true_type{};

template <typename T, typename enable = void>
struct is_row_vector_armadillo : std::false_type {};

template <typename T>
struct is_row_vector_armadillo<T,
	 ::pressio::mpl::enable_if_t<
	   std::is_same<T,
     	    arma::Row<typename T::elem_type>
			>::value
	   >
      > : std::true_type{};

template <typename T, typename enable = void>
struct is_vector_armadillo : std::false_type {};

template <typename T>
struct is_vector_armadillo<T,
	 ::pressio::mpl::enable_if_t<
	   is_row_vector_armadillo<T>::value or
	   is_column_vector_armadillo<T>::valu
	   >
      > : std::true_type{};

}}}//end namespace pressio::containers::meta
#endif
#endif
