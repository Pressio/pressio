
#ifndef ODE_IS_VALID_USER_DEFINED_OPS_IMPLICIT_ODE_HPP_
#define ODE_IS_VALID_USER_DEFINED_OPS_IMPLICIT_ODE_HPP_

#include "../../../core/src/meta/core_has_update_op_typedef.hpp"
#include "../../../core/src/meta/core_has_static_method_scale.hpp"
#include "../../../core/src/meta/core_has_static_method_add_to_diagonal.hpp"
#include "../../../core/src/meta/core_has_static_method_do_update_two_terms.hpp"

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
      ::rompp::core::meta::has_update_op_typedef<T>::value
      and
      ::rompp::core::meta::has_static_method_do_update_two_terms<
	typename T::update_op,
	scalar_t,
	typename core::details::traits<residual_t>::wrapped_t,
	typename core::details::traits<state_t>::wrapped_t,
	typename core::details::traits<state_t>::wrapped_t
	>::value
      and
      ::rompp::core::meta::has_static_method_scale<
	typename T::update_op,
	typename core::details::traits<jacobian_t>::wrapped_t,
	scalar_t
	>::value
      and
      ::rompp::core::meta::has_static_method_add_to_diagonal<
	typename T::update_op,
	typename core::details::traits<jacobian_t>::wrapped_t,
	scalar_t
	>::value
      >
  > : std::true_type{};
//--------------------------------------------


template<
  typename T,
  typename scalar_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t,
  typename enable = void
  >
struct is_valid_user_defined_ops_for_implicit_bdf2 : std::false_type{};

template<
  typename T,
  typename scalar_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t
  >
struct is_valid_user_defined_ops_for_implicit_bdf2<
  T, scalar_t, state_t, residual_t, jacobian_t,
    mpl::enable_if_t<
      ::rompp::core::meta::has_update_op_typedef<T>::value
      and
      ::rompp::core::meta::has_static_method_do_update_three_terms<
	typename T::update_op,
	scalar_t,
	typename core::details::traits<residual_t>::wrapped_t,
	typename core::details::traits<state_t>::wrapped_t,
	typename core::details::traits<state_t>::wrapped_t,
	typename core::details::traits<state_t>::wrapped_t
	>::value
      and
      ::rompp::core::meta::has_static_method_scale<
	typename T::update_op,
	typename core::details::traits<jacobian_t>::wrapped_t,
	scalar_t
	>::value
      and
      ::rompp::core::meta::has_static_method_add_to_diagonal<
	typename T::update_op,
	typename core::details::traits<jacobian_t>::wrapped_t,
	scalar_t
	>::value
      >
  > : std::true_type{};
//--------------------------------------------


template<
  typename T,
  typename scalar_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t,
  typename enable = void
  >
struct is_valid_user_defined_ops_for_implicit_ode : std::false_type{};

template<
  typename T,
  typename scalar_t,
  typename state_t,
  typename residual_t,
  typename jacobian_t
  >
struct is_valid_user_defined_ops_for_implicit_ode<
  T, scalar_t, state_t, residual_t, jacobian_t,
  mpl::enable_if_t<
    is_valid_user_defined_ops_for_implicit_euler<
      T, scalar_t, state_t, residual_t, jacobian_t
      >::value and
    is_valid_user_defined_ops_for_implicit_bdf2<
      T, scalar_t, state_t, residual_t, jacobian_t
      >::value
    >
  > : std::true_type{};

}}} // namespace rompp::ode::meta
#endif
