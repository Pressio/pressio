
#ifndef CORE_SHARED_TRAITS_HPP_
#define CORE_SHARED_TRAITS_HPP_

#include <type_traits>
#include "core_wrapped_types_enum.hpp"

namespace rompp{ namespace core{ namespace details {

//---------------------------------------
/// common traits of core containers
template<typename container_T,
	 typename wrapped_T,
	 bool is_vector_t,
	 bool is_matrix_t,
	 bool is_multi_vector_t,
	 WrappedPackageIdentifier wpid,
	 bool is_shared_mem_t>
struct containers_shared_traits{

  using wrapped_t = wrapped_T;
  using derived_t = container_T;

  static constexpr WrappedPackageIdentifier
  wrapped_package_identifier = wpid;

  static constexpr bool is_vector = is_vector_t;
  static constexpr bool is_matrix = is_matrix_t;
  static constexpr bool is_multi_vector = is_multi_vector_t;
  static constexpr bool is_shared_mem = is_shared_mem_t;

};


//---------------------------------------
/// common traits of core matrices
template<bool is_sparse_t>
struct matrix_shared_traits{

  static constexpr bool is_sparse = is_sparse_t;
  static constexpr bool is_dense = !is_sparse_t;

};


}}} // end namespace rompp::core::details
#endif
