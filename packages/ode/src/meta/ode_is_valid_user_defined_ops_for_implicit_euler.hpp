
#ifndef ODE_IS_VALID_USER_DEFINED_OPS_IMPLICIT_EULER_HPP_
#define ODE_IS_VALID_USER_DEFINED_OPS_IMPLICIT_EULER_HPP_

#include "../../../algebra/src/meta/algebra_has_update_op_typedef.hpp"
#include "../../../algebra/src/meta/algebra_has_static_method_scale.hpp"
#include "../../../algebra/src/meta/algebra_has_static_method_add_to_diagonal.hpp"
#include "../../../algebra/src/meta/algebra_has_static_method_do_update_two_terms.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<
  typename T,
  typename scalar_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t,
  typename enable = void
  >
struct is_valid_user_defined_ops_for_implicit_euler : std::false_type{};

template<
  typename T,
  typename scalar_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t
  >
struct is_valid_user_defined_ops_for_implicit_euler<
  T, scalar_t, state_t, residual_t, jacobian_t,
    mpl::enable_if_t<
      ::rompp::algebra::meta::has_update_op_typedef<T>::value
      and
      ::rompp::algebra::meta::has_static_method_do_update_two_terms<
	typename T::update_op,
	scalar_t,
	typename algebra::details::traits<residual_t>::wrapped_t,
	typename algebra::details::traits<state_t>::wrapped_t,
	typename algebra::details::traits<state_t>::wrapped_t
	>::value
      and
      ::rompp::algebra::meta::has_static_method_scale<
	typename T::update_op,
	typename algebra::details::traits<jacobian_t>::wrapped_t,
	scalar_t
	>::value
      and
      ::rompp::algebra::meta::has_static_method_add_to_diagonal<
	typename T::update_op,
	typename algebra::details::traits<jacobian_t>::wrapped_t,
	scalar_t
	>::value
      >
  > : std::true_type{};

}}} // namespace rompp::ode::meta
#endif
