
#ifndef SOLVERS_IS_LEGITIMATE_HESSIAN_FOR_GN_NORMEQ_HPP_
#define SOLVERS_IS_LEGITIMATE_HESSIAN_FOR_GN_NORMEQ_HPP_

#include "solvers_basic_meta.hpp"
#include "../base/solvers_linear_base.hpp"

namespace rompp{ namespace solvers{ namespace meta {

template <typename T, typename enable = void>
struct is_legitimate_hessian_for_gn_normeq
  : std::false_type{};

template <typename T>
struct is_legitimate_hessian_for_gn_normeq<
  T,
  ::rompp::mpl::enable_if_t<
    core::meta::is_core_matrix_wrapper<T>::value or
    core::meta::is_core_multi_vector_wrapper<T>::value
    >
  > : std::true_type{};

}}} // namespace rompp::solvers::meta
#endif
