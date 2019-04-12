
#ifndef SOLVERS_IS_LEGITIMATE_LINEAR_SOLVER_FOR_GN_NORMEQ_HPP_
#define SOLVERS_IS_LEGITIMATE_LINEAR_SOLVER_FOR_GN_NORMEQ_HPP_

#include "solvers_basic_meta.hpp"
#include "../base/solvers_linear_base.hpp"

namespace rompp{ namespace solvers{ namespace meta {

template <typename T, typename enable = void>
struct is_legitimate_linear_solver_for_gn_normeq
  : std::false_type{};

template <typename T>
struct is_legitimate_linear_solver_for_gn_normeq<
  T,
  ::rompp::mpl::enable_if_t<
    // the linear solver type has a public matrix_type typedef
    ::rompp::mpl::is_detected<has_matrix_typedef, T>::value and
    // the matrix_type is not void
    !std::is_void<typename T::matrix_type>::value and
    ::rompp::mpl::publicly_inherits_from<
      T,
      ::rompp::solvers::LinearBase<
	typename T::solver_t, typename T::matrix_type, T
	>
      >::value
    >
  > : std::true_type{};

}}} // namespace rompp::solvers::meta
#endif
