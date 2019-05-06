
#ifndef ODE_IS_VALID_USER_DEFINED_OPS_EXPLICIT_ODE_HPP_
#define ODE_IS_VALID_USER_DEFINED_OPS_EXPLICIT_ODE_HPP_

#include "../../../core/src/meta/core_has_update_op_typedef.hpp"
#include "../../../core/src/meta/core_has_static_method_do_update_one_term.hpp"
#include "../../../core/src/meta/core_has_static_method_do_update_four_terms.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t,
	 typename enable = void>
struct is_valid_user_defined_ops_for_explicit_euler : std::false_type{};

template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t>
struct is_valid_user_defined_ops_for_explicit_euler<
  T, scalar_t, state_t, residual_t,
    mpl::enable_if_t<
      ::rompp::core::meta::has_update_op_typedef<T>::value and
      ::rompp::core::meta::has_static_method_do_update_one_term<
	typename T::update_op,
	scalar_t,
	typename core::details::traits<state_t>::wrapped_t,
	typename core::details::traits<residual_t>::wrapped_t
	>::value
      >
  > : std::true_type{};
//--------------------------------------------


template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t,
	 typename enable = void>
struct is_valid_user_defined_ops_for_explicit_rk4 : std::false_type{};

template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t>
struct is_valid_user_defined_ops_for_explicit_rk4<
  T, scalar_t, state_t, residual_t,
    mpl::enable_if_t<
      ::rompp::core::meta::has_update_op_typedef<T>::value and
      ::rompp::core::meta::has_static_method_do_update_four_terms<
    typename T::update_op,
	scalar_t,
	typename core::details::traits<state_t>::wrapped_t,
	typename core::details::traits<residual_t>::wrapped_t,
	typename core::details::traits<residual_t>::wrapped_t,
	typename core::details::traits<residual_t>::wrapped_t,
	typename core::details::traits<residual_t>::wrapped_t
	>::value
      >
  > : std::true_type{};
//--------------------------------------------


template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t,
	 typename enable = void>
struct is_valid_user_defined_ops_for_explicit_ode : std::false_type{};

template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t>
struct is_valid_user_defined_ops_for_explicit_ode<
  T, scalar_t, state_t, residual_t,
    mpl::enable_if_t<
      is_valid_user_defined_ops_for_explicit_euler<
	T, scalar_t, state_t, residual_t
      	>::value and
      is_valid_user_defined_ops_for_explicit_rk4<
	T, scalar_t, state_t, residual_t
      	>::value
      >
  > : std::true_type{};

}}} // namespace rompp::ode::meta
#endif
