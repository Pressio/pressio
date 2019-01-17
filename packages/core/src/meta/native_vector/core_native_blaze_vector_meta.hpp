
#ifdef HAVE_BLAZE
#ifndef CORE_NATIVE_BLAZE_VECTOR_META_HPP_
#define CORE_NATIVE_BLAZE_VECTOR_META_HPP_

#include "../core_meta_basic.hpp"
#include <blaze/math/DynamicVector.h>
#include <blaze/math/StaticVector.h>


namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_static_vector_blaze : std::false_type {};

template <typename T>
struct is_static_vector_blaze<T,
	 core::meta::enable_if_t<
	   blaze::IsStatic<typename T::This>::value
	   >
      > : std::true_type{};

template <typename T, typename enable = void>
struct is_dynamic_row_vector_blaze : std::false_type {};

template <typename T>
struct is_dynamic_row_vector_blaze<T,
	 core::meta::enable_if_t<
	   std::is_same<T, blaze::DynamicVector<
				typename T::ElementType,
				blaze::rowVector>
			>::value
	   >
      > : std::true_type{};

template <typename T, typename enable = void>
struct is_dynamic_column_vector_blaze : std::false_type {};

template <typename T>
struct is_dynamic_column_vector_blaze<T,
	 core::meta::enable_if_t<
	   std::is_same<T, blaze::DynamicVector<
				typename T::ElementType,
				blaze::columnVector>
			>::value
	   >
      > : std::true_type{};

template <typename T, typename enable = void>
struct is_dynamic_vector_blaze : std::false_type {};

template <typename T>
struct is_dynamic_vector_blaze<T,
	   core::meta::enable_if_t<
	     is_dynamic_row_vector_blaze<T>::value ||
	     is_dynamic_column_vector_blaze<T>::value
	   >
      > : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
#endif
