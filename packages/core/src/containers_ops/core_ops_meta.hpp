
#ifndef CORE_OPS_META_HPP_
#define CORE_OPS_META_HPP_

#include "../matrix/core_matrix_meta.hpp"
#include "../vector/core_vector_meta.hpp"
#include "../multi_vector/core_multi_vector_meta.hpp"


namespace rompp{
namespace core{
namespace meta {

template <typename T, typename T2, typename enable = void>
struct wrappers_have_same_scalar : std::false_type {};

template <typename T, typename T2>
struct wrappers_have_same_scalar<T, T2,
  core::meta::enable_if_t<
    std::is_same<typename T::scalar_t,
		 typename T2::scalar_t
		 >::value
    >
  > : std::true_type{};
  
 
} // namespace meta
} // namespace core

}//end namespace rompp
#endif
