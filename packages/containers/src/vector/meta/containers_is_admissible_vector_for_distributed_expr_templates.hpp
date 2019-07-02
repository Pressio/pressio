
#ifndef CONTAINERS_IS_ADMISSIBLE_VECTOR_FOR_DISTRIBUTED_EXP_TEMPL_HPP_
#define CONTAINERS_IS_ADMISSIBLE_VECTOR_FOR_DISTRIBUTED_EXP_TEMPL_HPP_

#include "containers_is_vector_wrapper.hpp"

namespace pressio{ namespace containers{ namespace meta {

template <typename T,
    typename enable = void>
struct is_admissible_vec_for_dist_expression : std::false_type{};

template <typename T>
struct is_admissible_vec_for_dist_expression<T,
      ::pressio::mpl::enable_if_t<
  containers::meta::is_vector_wrapper<T>::value &&
  !containers::details::traits<T>::is_shared_mem
      >> : std::true_type{};

}}}//end namespace pressio::containers::meta
#endif
