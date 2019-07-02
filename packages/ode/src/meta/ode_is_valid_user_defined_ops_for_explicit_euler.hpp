
#ifndef ODE_IS_VALID_USER_DEFINED_OPS_EXPLICIT_EULER_HPP_
#define ODE_IS_VALID_USER_DEFINED_OPS_EXPLICIT_EULER_HPP_

#include "../../../containers/src/meta/containers_has_update_op_typedef.hpp"
#include "../../../containers/src/meta/containers_has_static_method_do_update_one_term.hpp"

namespace pressio{ namespace ode{ namespace meta {

template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t,
	 typename enable = void>
struct is_valid_user_defined_ops_for_explicit_euler
  : std::false_type{};


template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t>
struct is_valid_user_defined_ops_for_explicit_euler<
  T, scalar_t, state_t, residual_t,
    mpl::enable_if_t<
      ::pressio::containers::meta::has_update_op_typedef<T>::value and
      ::pressio::containers::meta::is_vector_wrapper<state_t>::value and
      ::pressio::containers::meta::has_static_method_do_update_one_term<
	typename T::update_op,
	scalar_t,
	typename containers::details::traits<state_t>::wrapped_t,
	typename containers::details::traits<residual_t>::wrapped_t
	>::value
      >
  > : std::true_type{};


#ifdef HAVE_PYBIND11
template<typename T,
	 typename scalar_t,
	 typename state_t,
	 typename residual_t>
struct is_valid_user_defined_ops_for_explicit_euler<
  T, scalar_t, state_t, residual_t,
    mpl::enable_if_t<
      ::pressio::containers::meta::is_cstyle_array_pybind11<state_t>::value and
      ::pressio::containers::meta::is_cstyle_array_pybind11<residual_t>::value and
      ::pressio::containers::meta::has_update_op_typedef<T>::value and
      ::pressio::containers::meta::has_static_method_do_update_one_term<
	typename T::update_op,
	scalar_t, state_t, residual_t
	>::value
      >
  > : std::true_type{};
#endif

}}} // namespace pressio::ode::meta
#endif
