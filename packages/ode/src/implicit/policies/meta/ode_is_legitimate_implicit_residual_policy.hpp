
#ifndef ODE_POLICIES_META_IS_LEGITIMATE_IMPLICIT_RESIDUAL_POLICIES_HPP_
#define ODE_POLICIES_META_IS_LEGITIMATE_IMPLICIT_RESIDUAL_POLICIES_HPP_

#include "../base/ode_implicit_residual_policy_base.hpp"

namespace rompp{ namespace ode{ namespace meta {

// template<typename policy_t, typename enable = void>
// struct inherits_from_implicit_residual_pol_base : std::false_type{};

// template <typename policy_t>
// struct inherits_from_implicit_residual_pol_base<
//   policy_t,
//   typename std::enable_if<
//     core::meta::publicly_inherits_from<
//       policy_t,
//       ode::policy::ImplicitResidualPolicyBase<policy_t>
//       >::value
//     >::type
//   > : std::true_type{};
// //------------------------------------------------------------------


template <typename T, ::rompp::ode::ImplicitEnum odeMethod,
	  int numAuxStates, typename t1, typename t2, typename t3,
	  typename enable = void>
struct is_residual_callable_with_five_args : std::false_type
{};

template <typename T, ::rompp::ode::ImplicitEnum odeMethod,
	  int numAuxStates, typename t1, typename t2, typename t3>
struct is_residual_callable_with_five_args<
  T, odeMethod, numAuxStates, t1, t2, t3,
  core::meta::enable_if_t<
    !std::is_void<
      decltype(std::declval<T>().template operator()<
	   odeMethod,
	   numAuxStates>( std::declval<const t1 &>(),
			  std::declval<const std::array<t1,numAuxStates> &>(),
			  std::declval<const t2&>(),
			  std::declval<t3>(),
			  std::declval<t3>()
			  )
	       )
      >::value
    >
  > : std::true_type
{};
//--------------------------------------------------------------------

template <typename T, ::rompp::ode::ImplicitEnum odeMethod,
	  int numAuxStates, typename t1, typename t2,
	  typename t3, typename t4, typename enable = void>
struct is_residual_callable_with_six_args : std::false_type
{};

template <typename T, ::rompp::ode::ImplicitEnum odeMethod,
	  int numAuxStates, typename t1, typename t2,
	  typename t3, typename t4>
struct is_residual_callable_with_six_args<
  T, odeMethod, numAuxStates, t1, t2, t3, t4,
  core::meta::enable_if_t<
    std::is_void<
      decltype(std::declval<T>().template operator()<
	       odeMethod,
	       numAuxStates>( std::declval<const t1 &>(),
			std::declval<t2 &>(),
			std::declval<const std::array<t1,numAuxStates> &>(),
			std::declval<const t3&>(),
			std::declval<t4>(),
			std::declval<t4>()
			)
	       )
      >::value
    >
  > : std::true_type
{};
//-------------------------------------------------------------------

template<ImplicitEnum name, int nstates,
	 typename T, typename state_t, typename residual_t,
	 typename model_t, typename scalar_t, typename enable = void>
struct is_legitimate_implicit_residual_policy
  : std::false_type{};

template<ImplicitEnum name, int nstates,
	 typename T, typename state_t, typename residual_t,
	 typename model_t, typename scalar_t>
struct is_legitimate_implicit_residual_policy
<name, nstates, T, state_t, residual_t, model_t, scalar_t,
 core::meta::enable_if_t<
   // inherits from base
   /*inherits_from_implicit_residual_pol_base<T>::value and*/
   // is callable with 6 args
   is_residual_callable_with_six_args<
     T, name, nstates, state_t, residual_t, model_t, scalar_t
     >::value and
   // is callable with 5 args
   is_residual_callable_with_five_args<
     T, name, nstates, state_t, model_t, scalar_t
     >::value
   >
 > : std::true_type{};
//------------------------------------------------------------------

template<typename T, typename ... args>
using is_legitimate_implicit_euler_residual_policy =
  is_legitimate_implicit_residual_policy<
  ImplicitEnum::Euler, 1, T, args...>;

template<typename T, typename ... args>
using is_legitimate_implicit_bdf2_residual_policy =
  is_legitimate_implicit_residual_policy<
  ImplicitEnum::BDF2, 2, T, args...>;


}}} // namespace rompp::ode::meta
#endif
