
#ifndef ODE_IS_VALID_USER_DEFINED_OPS_EXPLICIT_ODE_HPP_
#define ODE_IS_VALID_USER_DEFINED_OPS_EXPLICIT_ODE_HPP_

#include "ode_basic_meta.hpp"

namespace rompp{ namespace ode{ namespace meta {

/*
  User-defined ops for explicit ode are valid if:
  * it contains a typedef for updating a container
  * the above typedef has a do_update static method
 */

template <typename T, typename = void>
struct has_update_op_typedef : std::false_type{};

template <typename T>
struct has_update_op_typedef<
  T, mpl::enable_if_t<
       !std::is_void<
	 typename T::update_op
	 >::value
       >
  > : std::true_type{};
//--------------------------------------------


template <typename T,
	  typename scalar_t,
	  typename state_t,
	  typename residual_t,
	  typename = void>
struct has_static_do_update_one_term : std::false_type{};

template <typename T,
	  typename sc_t,
	  typename state_t,
	  typename residual_t >
struct has_static_do_update_one_term<
  T, sc_t, state_t, residual_t,
  mpl::enable_if_t<
    std::is_void<
      decltype
      (
       T::do_update
       (
	std::declval<typename core::details::traits<state_t>::wrapped_t &>(),
	std::declval<const sc_t>(),
	std::declval<const typename core::details::traits<residual_t>::wrapped_t &>(),
	std::declval<const sc_t>()
	)
       )
      >::value
    >
  > : std::true_type{};
//--------------------------------------------


// template <typename T,
// 	  typename scalar_t,
// 	  typename state_t,
// 	  typename residual_t,
// 	  typename = void>
// struct has_static_do_update_four_terms : std::false_type{};

// template <typename T,
// 	  typename sc_t,
// 	  typename state_t,
// 	  typename residual_t >
// struct has_static_do_update_four_terms<
//   T, sc_t, state_t, residual_t,
//   mpl::enable_if_t<
//     mpl::void_t<
//       decltype(
// 	       T::do_update(
// 			    std::declval<state_t &>(),
// 			    std::declval<const sc_t >(),
// 			    std::declval<residual_t &>(),
// 			    std::declval<const sc_t >(),
// 			    std::declval<residual_t &>(),
// 			    std::declval<const sc_t >(),
// 			    std::declval<residual_t &>(),
// 			    std::declval<const sc_t >(),
// 			    std::declval<residual_t &>(),
// 			    std::declval<const sc_t >()
// 			    )
// 		)
//       >
//     >
//   > : std::true_type{};
// //--------------------------------------------


template <typename T,
	  typename scalar_t,
	  typename state_t,
	  typename residual_t,
	  typename = void>
struct update_op_has_all_needed_methods : std::false_type{};

template <typename T,
	  typename scalar_t,
	  typename state_t,
	  typename residual_t >
struct update_op_has_all_needed_methods<
  T, scalar_t, state_t, residual_t,
  mpl::enable_if_t<
    has_static_do_update_one_term<T, scalar_t, state_t, residual_t>::value
    //has_static_do_update_four_terms<T, scalar_t, state_t, residual_t>::value
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
      ::rompp::ode::meta::has_update_op_typedef<T>::value and
      ::rompp::ode::meta::update_op_has_all_needed_methods<
      	typename T::update_op, scalar_t, state_t, residual_t
      	>::value
      >
  > : std::true_type{};

}}} // namespace rompp::ode::meta
#endif
