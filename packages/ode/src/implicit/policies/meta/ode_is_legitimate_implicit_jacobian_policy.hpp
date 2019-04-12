
#ifndef ODE_POLICIES_META_IS_LEGITIMATE_IMPLICIT_JACOBIAN_POLICIES_HPP_
#define ODE_POLICIES_META_IS_LEGITIMATE_IMPLICIT_JACOBIAN_POLICIES_HPP_

#include "../base/ode_jacobian_policy_base.hpp"

namespace rompp{ namespace ode{ namespace meta {

// template<typename policy_t, typename enable = void>
// struct inherits_from_jacobian_pol_base : std::false_type{};

// template <typename policy_t>
// struct inherits_from_jacobian_pol_base<
//   policy_t,
//   typename std::enable_if<
//     core::meta::publicly_inherits_from<
//       policy_t,
//       ode::policy::JacobianPolicyBase<policy_t>
//       >::value
//     >::type
//   > : std::true_type{};
// //------------------------------------------------------------------


template <typename T, ::rompp::ode::ImplicitEnum odeMethod,
	  typename t1, typename t2, typename t3,
	  typename enable = void>
struct is_jacobian_callable_with_four_args : std::false_type
{};

template <typename T, ::rompp::ode::ImplicitEnum odeMethod,
	  typename t1, typename t2, typename t3>
struct is_jacobian_callable_with_four_args<
  T, odeMethod, t1, t2, t3,
  ::rompp::mpl::enable_if_t<
    !std::is_void<
      decltype(std::declval<T>().template operator()<
	       odeMethod>( std::declval<const t1 &>(),
			  std::declval<const t2&>(),
			  std::declval<t3>(),
			  std::declval<t3>()
			  )
	       )
      >::value
    >
  > : std::true_type
{};
//---------------------------------------------------------------------

template <typename T, ::rompp::ode::ImplicitEnum odeMethod,
	  typename t1, typename t2, typename t3, typename t4,
	  typename enable = void>
struct is_jacobian_callable_with_five_args : std::false_type
{};

template <typename T, ::rompp::ode::ImplicitEnum odeMethod,
	  typename t1, typename t2, typename t3, typename t4>
struct is_jacobian_callable_with_five_args<
  T, odeMethod, t1, t2, t3, t4,
  ::rompp::mpl::enable_if_t<
    std::is_void<
      decltype(std::declval<T>().template operator()<
	       odeMethod>( std::declval<const t1 &>(),
			   std::declval<t2 &>(),
			   std::declval<const t3&>(),
			   std::declval<t4>(),
			   std::declval<t4>()
			   )
	       )
      >::value
    >
  > : std::true_type
{};
//---------------------------------------------------------------------


template<ImplicitEnum name, typename T,
	 typename state_t, typename jacobian_t,
	 typename model_t, typename scalar_t, typename enable = void>
struct is_legitimate_implicit_jacobian_policy : std::false_type
{};

template<ImplicitEnum name, typename T,
	 typename state_t, typename jacobian_t,
	 typename model_t, typename scalar_t>
struct is_legitimate_implicit_jacobian_policy
<name, T, state_t, jacobian_t, model_t, scalar_t,
 ::rompp::mpl::enable_if_t<
   /*inherits_from_jacobian_pol_base<T>::value and*/
   is_jacobian_callable_with_four_args<
     T, name, state_t, model_t, scalar_t>::value and
   is_jacobian_callable_with_five_args<
     T, name, state_t, jacobian_t,
     model_t, scalar_t>::value
   >
 > : std::true_type{};
//------------------------------------------------------------------

template<typename T, typename ... args>
using is_legitimate_implicit_euler_jacobian_policy =
  is_legitimate_implicit_jacobian_policy<
  ::rompp::ode::ImplicitEnum::Euler, T, args...>;

template<typename T, typename ... args>
using is_legitimate_implicit_bdf2_jacobian_policy =
  is_legitimate_implicit_jacobian_policy<
  ::rompp::ode::ImplicitEnum::BDF2, T, args...>;


}}} // namespace rompp::ode::meta
#endif
