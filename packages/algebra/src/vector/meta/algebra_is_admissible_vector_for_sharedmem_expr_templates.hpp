
#ifndef ALGEBRA_IS_ADMISSIBLE_VECTOR_FOR_SHAREDMEM_EXP_TEMPL_HPP_
#define ALGEBRA_IS_ADMISSIBLE_VECTOR_FOR_SHAREDMEM_EXP_TEMPL_HPP_

#include "algebra_is_algebra_vector_wrapper.hpp"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T,
    typename enable = void>
struct is_admissible_vec_for_sharedmem_expression : std::false_type{};

template <typename T>
struct is_admissible_vec_for_sharedmem_expression<T,
      ::rompp::mpl::enable_if_t<
  algebra::meta::is_algebra_vector_wrapper<T>::value &&
  algebra::details::traits<T>::is_shared_mem
      >> : std::true_type{};

}}}//end namespace rompp::algebra::meta
#endif
