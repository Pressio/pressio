
#ifndef SOLVERS_IS_LEGITIMATE_QR_SOLVER_FOR_GN_QR_HPP_
#define SOLVERS_IS_LEGITIMATE_QR_SOLVER_FOR_GN_QR_HPP_

#include "solvers_basic_meta.hpp"
#include "../base/solvers_linear_base.hpp"
#include "../../../qr/src/qr_fwd.hpp"
#include "../../../qr/src/qr_traits.hpp"

namespace pressio{ namespace solvers{ namespace meta {

template <typename T, typename enable = void>
struct is_legitimate_qr_solver_for_gn_qr
  : std::false_type{};

template <typename T>
struct is_legitimate_qr_solver_for_gn_qr<
  T,
  ::pressio::mpl::void_t<
    typename ::pressio::qr::details::traits<T>::concrete_t
    >
  > : std::true_type{};

}}} // namespace pressio::solvers::meta
#endif
