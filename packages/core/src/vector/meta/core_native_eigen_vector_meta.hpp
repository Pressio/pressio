
#ifndef CORE_NATIVE_EIGEN_VECTOR_META_HPP_
#define CORE_NATIVE_EIGEN_VECTOR_META_HPP_

#include "../../meta/core_meta_basic.hpp"
#include <Eigen/Dense>

namespace rompp{ namespace core{ namespace meta {


template <typename T, typename enable = void>
struct is_dynamic_row_vector_eigen : std::false_type {};

template <typename T>
struct is_dynamic_row_vector_eigen<T,
	 ::rompp::mpl::enable_if_t<
	   std::is_same<T,
			Eigen::Matrix<typename T::Scalar,
				      1, Eigen::Dynamic>
			>::value
	   >
      > : std::true_type{};
//----------------------------------------------

template <typename T, typename enable = void>
struct is_static_row_vector_eigen : std::false_type {};

template <typename T>
struct is_static_row_vector_eigen<T,
	 ::rompp::mpl::enable_if_t<
	   std::is_same<T,
			Eigen::Matrix<typename T::Scalar,
				      1,
				      T::ColsAtCompileTime>
			>::value and
	   !is_dynamic_row_vector_eigen<T>::value
	   >
      > : std::true_type{};
//----------------------------------------------


template <typename T, typename enable = void>
struct is_dynamic_column_vector_eigen : std::false_type {};

template <typename T>
struct is_dynamic_column_vector_eigen<T,
	 ::rompp::mpl::enable_if_t<
	   std::is_same<T,
			Eigen::Matrix<typename T::Scalar,
				      Eigen::Dynamic, 1>
			>::value
	   >
      > : std::true_type{};
//----------------------------------------------

template <typename T, typename enable = void>
struct is_static_column_vector_eigen : std::false_type {};

template <typename T>
struct is_static_column_vector_eigen<T,
	 ::rompp::mpl::enable_if_t<
	   std::is_same<T,
			Eigen::Matrix<typename T::Scalar,
				      T::RowsAtCompileTime,
				      1>
			>::value and
	   !is_dynamic_column_vector_eigen<T>::value
	   >
      > : std::true_type{};
//----------------------------------------------


template <typename T, typename enable = void>
struct is_static_vector_eigen : std::false_type {};

template <typename T>
struct is_static_vector_eigen<T,
	 ::rompp::mpl::enable_if_t<
	   is_static_row_vector_eigen<T>::value ||
	   is_static_column_vector_eigen<T>::value
	   >
      > : std::true_type{};
//----------------------------------------------


template <typename T, typename enable = void>
struct is_dynamic_vector_eigen : std::false_type {};

template <typename T>
struct is_dynamic_vector_eigen<T,
	 ::rompp::mpl::enable_if_t<
	   is_dynamic_row_vector_eigen<T>::value or
	   is_dynamic_column_vector_eigen<T>::value
	   >
      > : std::true_type{};
//----------------------------------------------


template <typename T, typename enable = void>
struct is_vector_eigen : std::false_type {};

template <typename T>
struct is_vector_eigen< T,
	 ::rompp::mpl::enable_if_t<
	   is_dynamic_vector_eigen<T>::value or
	   is_static_vector_eigen<T>::value
	   >
     > : std::true_type{};


}}}//end namespace rompp::core::meta
#endif
