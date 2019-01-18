
#ifndef SOLVERS_IMPL_HESSIAN_APPROX_HELPERS_MIXIN_HPP
#define SOLVERS_IMPL_HESSIAN_APPROX_HELPERS_MIXIN_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../../../CORE_OPS"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template<typename J_t, typename enable = void>
struct HessianApproxHelper;


// when J is matrix wrapper, the hessian J^T*J
// is computed by doing product of J^T*J
template<typename J_t>
struct HessianApproxHelper<
  J_t,
  core::meta::enable_if_t<
    core::meta::is_core_matrix_wrapper<J_t>::value
    >>
{
  template <typename result_t>
  void operator()(J_t & J, result_t & result) const{
    constexpr bool transposeJ = true;
    ::rompp::core::ops::product<J_t, J_t, result_t,
				transposeJ>(J, J, result);
  }

  auto operator()(J_t & J) const
    -> decltype(::rompp::core::ops::product<J_t, J_t, true>(J, J)) {
    return ::rompp::core::ops::product<J_t, J_t, true>(J, J);
  }
};


// when J is multivector wrapper, the hessian J^T*J
// is computed by doing the DOT of J*J
template<typename J_t>
struct HessianApproxHelper<
  J_t,
  core::meta::enable_if_t<
    core::meta::is_core_multi_vector_wrapper<J_t>::value
    >>
{
  template <typename result_t>
  void operator()(J_t & J, result_t & result) const{
    constexpr bool transposeJ = true;
    ::rompp::core::ops::dot(J, J, result);
  }
};


}}}} //end namespace rompp::solvers::iterative::impl
#endif
