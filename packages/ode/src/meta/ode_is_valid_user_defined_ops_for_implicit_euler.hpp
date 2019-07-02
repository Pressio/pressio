
#ifndef ODE_IS_VALID_USER_DEFINED_OPS_IMPLICIT_EULER_HPP_
#define ODE_IS_VALID_USER_DEFINED_OPS_IMPLICIT_EULER_HPP_

#include "../../../containers/src/meta/containers_has_update_op_typedef.hpp"
#include "../../../containers/src/meta/containers_has_static_method_scale.hpp"
#include "../../../containers/src/meta/containers_has_static_method_add_to_diagonal.hpp"
#include "../../../containers/src/meta/containers_has_static_method_do_update_two_terms.hpp"

namespace pressio{ namespace ode{ namespace meta {

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
      ::pressio::containers::meta::has_update_op_typedef<T>::value
      and
      ::pressio::containers::meta::has_static_method_do_update_two_terms<
	typename T::update_op,
	scalar_t,
	typename containers::details::traits<residual_t>::wrapped_t,
	typename containers::details::traits<state_t>::wrapped_t,
	typename containers::details::traits<state_t>::wrapped_t
	>::value
      and
      ::pressio::containers::meta::has_static_method_scale<
	typename T::update_op,
	typename containers::details::traits<jacobian_t>::wrapped_t,
	scalar_t
	>::value
      and
      ::pressio::containers::meta::has_static_method_add_to_diagonal<
	typename T::update_op,
	typename containers::details::traits<jacobian_t>::wrapped_t,
	scalar_t
	>::value
      >
  > : std::true_type{};

}}} // namespace pressio::ode::meta
#endif
