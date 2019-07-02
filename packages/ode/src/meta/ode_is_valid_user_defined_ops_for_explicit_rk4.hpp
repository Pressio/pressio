
#ifndef ODE_IS_VALID_USER_DEFINED_OPS_EXPLICIT_RK4_HPP_
#define ODE_IS_VALID_USER_DEFINED_OPS_EXPLICIT_RK4_HPP_

#include "../../../containers/src/meta/containers_has_update_op_typedef.hpp"
#include "../../../containers/src/meta/containers_has_static_method_do_update_two_terms.hpp"
#include "../../../containers/src/meta/containers_has_static_method_do_update_four_terms.hpp"

namespace pressio{ namespace ode{ namespace meta {

template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t,
	 typename enable = void>
struct is_valid_user_defined_ops_for_explicit_rk4
  : std::false_type{};

template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t>
struct is_valid_user_defined_ops_for_explicit_rk4<
  T, scalar_t, state_t, residual_t,
    mpl::enable_if_t<
      ::pressio::containers::meta::is_vector_wrapper<state_t>::value
      and
      ::pressio::containers::meta::has_update_op_typedef<T>::value
      and
      ::pressio::containers::meta::has_static_method_do_update_two_terms<
	typename T::update_op,
	scalar_t,
	typename containers::details::traits<state_t>::wrapped_t,
	typename containers::details::traits<residual_t>::wrapped_t,
	typename containers::details::traits<residual_t>::wrapped_t
	>::value
      and
      ::pressio::containers::meta::has_static_method_do_update_four_terms<
	typename T::update_op,
	scalar_t,
	typename containers::details::traits<state_t>::wrapped_t,
	typename containers::details::traits<residual_t>::wrapped_t,
	typename containers::details::traits<residual_t>::wrapped_t,
	typename containers::details::traits<residual_t>::wrapped_t,
	typename containers::details::traits<residual_t>::wrapped_t
	>::value
      >
  > : std::true_type{};


#ifdef HAVE_PYBIND11
template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t>
struct is_valid_user_defined_ops_for_explicit_rk4<
  T, scalar_t, state_t, residual_t,
    mpl::enable_if_t<
      ::pressio::containers::meta::is_array_pybind11<state_t>::value
      and
      ::pressio::containers::meta::is_array_pybind11<residual_t>::value
      and
      ::pressio::containers::meta::has_update_op_typedef<T>::value
      and
      ::pressio::containers::meta::has_static_method_do_update_two_terms<
	typename T::update_op,
	scalar_t, state_t, residual_t, residual_t
	>::value
      and
      ::pressio::containers::meta::has_static_method_do_update_four_terms<
	typename T::update_op,
	scalar_t, state_t, residual_t, residual_t, residual_t, residual_t
	>::value
      >
  > : std::true_type{};
#endif


}}} // namespace pressio::ode::meta
#endif
