
#ifndef SOLVERS_IMPL_HESSIAN_APPROX_HELPERS_POLICY_HPP
#define SOLVERS_IMPL_HESSIAN_APPROX_HELPERS_POLICY_HPP

#include "../../solvers_ConfigDefs.hpp"
#include "../../../../ALGEBRA_OPS"

namespace rompp{ namespace solvers{ namespace iterative{ namespace impl{

template<typename J_t, typename enable = void>
struct HessianApproxHelper;


// when J is matrix wrapper, the hessian J^T*J
// is computed by doing product of J^T*J
template<typename J_t>
struct HessianApproxHelper<
  J_t,
  ::rompp::mpl::enable_if_t<
    algebra::meta::is_algebra_matrix_wrapper<J_t>::value
    >>{

  template <typename result_t>
  static void evaluate(J_t & J, result_t & result)
  {
    constexpr bool transposeJ = true;
    ::rompp::algebra::ops::product<J_t, J_t, result_t,
				transposeJ>(J, J, result);
  }

  static auto evaluate(J_t & J)
    -> decltype(::rompp::algebra::ops::product<J_t, J_t, true>(J, J))
  {
    return ::rompp::algebra::ops::product<J_t, J_t, true>(J, J);
  }
};


/* when J is multivector wrapper, the hessian J^T*J
 * is computed by doing the J self_dot J
 * this ensures that we leverage symmetry of the result,
 * since self_dot computes only half of the result matrix
 * and fills the rest by symmetry
 */
template<typename J_t>
struct HessianApproxHelper<
  J_t,
  ::rompp::mpl::enable_if_t<
    algebra::meta::is_algebra_multi_vector_wrapper<J_t>::value
    >>
{
  template <typename result_t>
  static void evaluate(J_t & J, result_t & result) {
    ::rompp::algebra::ops::dot_self(J, result);
  }

  static auto evaluate(J_t & J)
    -> decltype(::rompp::algebra::ops::dot_self(J)) {
    return ::rompp::algebra::ops::dot_self(J);
  }
};


}}}} //end namespace rompp::solvers::iterative::impl
#endif
