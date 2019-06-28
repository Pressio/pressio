
#ifndef ALGEBRA_SHARED_TRAITS_HPP_
#define ALGEBRA_SHARED_TRAITS_HPP_

#include "algebra_wrapped_types_enum.hpp"

namespace rompp{ namespace algebra{ namespace details {

//---------------------------------------
/// common traits of algebra containers
template<typename container_T,
	 typename wrapped_T,
	 bool is_vector_t,
	 bool is_matrix_t,
	 bool is_multi_vector_t,
	 WrappedPackageIdentifier wpid,
	 bool is_shared_mem_t,
	 bool is_static_t>
struct containers_shared_traits{

  using wrapped_t = wrapped_T;
  using derived_t = container_T;

  static constexpr WrappedPackageIdentifier
  wrapped_package_identifier = wpid;

  static constexpr bool is_vector = is_vector_t;
  static constexpr bool is_matrix = is_matrix_t;
  static constexpr bool is_multi_vector = is_multi_vector_t;
  static constexpr bool is_shared_mem = is_shared_mem_t;
  static constexpr bool is_static = is_static_t;
  static constexpr bool is_dynamic = !is_static_t;

  // by default, any container is not admissible to expr templates
  // the ones that are, will overwrite this
  static constexpr bool is_admissible_for_expression_templates = false;
};


//---------------------------------------
/// common traits of algebra matrices
template<bool is_sparse_t>
struct matrix_shared_traits{

  static constexpr bool is_sparse = is_sparse_t;
  static constexpr bool is_dense = !is_sparse_t;

};


}}} // end namespace rompp::algebra::details
#endif
