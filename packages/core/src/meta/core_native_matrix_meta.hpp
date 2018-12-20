
#ifndef CORE_NATIVE_MATRIX_MATRIX_META_HPP_
#define CORE_NATIVE_MATRIX_MATRIX_META_HPP_

#include "core_meta_basic.hpp"
// #include "core_native_vector_meta.hpp"
#include "core_native_trilinos_matrix_meta.hpp"
#include "core_native_eigen_matrix_meta.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_matrix_dense_sharedmem_stdlib : std::false_type {};

template <typename T>
struct is_matrix_dense_sharedmem_stdlib<T,
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
