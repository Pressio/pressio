
#ifndef CORE_IS_ADMISSIBLE_VECTOR_FOR_DISTRIBUTED_EXP_TEMPL_HPP_
#define CORE_IS_ADMISSIBLE_VECTOR_FOR_DISTRIBUTED_EXP_TEMPL_HPP_

#include "core_is_core_vector_wrapper.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T,
    typename enable = void>
struct is_admissible_vec_for_dist_expression : std::false_type{};

template <typename T>
struct is_admissible_vec_for_dist_expression<T,
      ::rompp::mpl::enable_if_t<
  core::meta::is_core_vector_wrapper<T>::value &&
  !core::details::traits<T>::is_shared_mem
      >> : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
