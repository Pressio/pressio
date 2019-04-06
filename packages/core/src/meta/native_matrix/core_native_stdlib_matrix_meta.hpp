
#ifndef CORE_NATIVE_STDLIB_MATRIX_META_HPP_
#define CORE_NATIVE_STDLIB_MATRIX_META_HPP_

#include "../core_meta_basic.hpp"
#include <vector>

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_dense_matrix_stdlib : std::false_type {};

template <typename T>
struct is_dense_matrix_stdlib<T,
     typename
     std::enable_if<
       std::is_same<T,std::vector<
			std::vector<typename
			T::value_type::value_type
			>
		     >
		    >::value
		    >::type
    > : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
