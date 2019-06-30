
#ifndef CONTAINERS_IS_ADMISSIBLE_VECTOR_FOR_SHAREDMEM_EXP_TEMPL_HPP_
#define CONTAINERS_IS_ADMISSIBLE_VECTOR_FOR_SHAREDMEM_EXP_TEMPL_HPP_

#include "containers_is_vector_wrapper.hpp"

namespace rompp{ namespace containers{ namespace meta {

template <typename T,
    typename enable = void>
struct is_admissible_vec_for_sharedmem_expression : std::false_type{};

template <typename T>
struct is_admissible_vec_for_sharedmem_expression<T,
      ::rompp::mpl::enable_if_t<
  containers::meta::is_vector_wrapper<T>::value &&
  containers::details::traits<T>::is_shared_mem
      >> : std::true_type{};

}}}//end namespace rompp::containers::meta
#endif
