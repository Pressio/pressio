

#ifndef ODE_ode_admissible_implicit_jacobian_policy_regular_stepper_HPP_
#define ODE_ode_admissible_implicit_jacobian_policy_regular_stepper_HPP_

namespace pressio{ namespace ode{ namespace meta {

template<
  typename T,
  typename tag,
  typename state_t,
  typename jacobian_t,
  typename system_t,
  typename scalar_t,
  typename enable = void
  >
struct admissible_implicit_jacobian_policy_regular_stepper : std::false_type
{};


template<
  typename T,
  typename tag,
  typename state_t,
  typename jacobian_t,
  typename system_t,
  typename scalar_t
  >
struct admissible_implicit_jacobian_policy_regular_stepper
<T, tag, state_t, jacobian_t, system_t, scalar_t,
 ::pressio::mpl::enable_if_t<
   std::is_same<
     jacobian_t,
     decltype
     (
      std::declval<T const>().operator()
      (
       std::declval<state_t const &>(),
       std::declval<system_t const &>()
       )
      )
     >::value
   and
   std::is_void<
     decltype
     (
      std::declval<T const>().template operator()
      <tag>(
	    std::declval<state_t const &>(),
	    std::declval<system_t const &>(),
	    std::declval<scalar_t const &>(),
	    std::declval<scalar_t const &>(),
	    std::declval<::pressio::ode::types::step_t const &>(),
	    std::declval<jacobian_t &>()
	)
      )
   >::value
   >
 > : std::true_type{};
//------------------------------------------------------------------

template<typename T, typename ... args>
using admissible_implicit_euler_jacobian_policy_regular_stepper =
  admissible_implicit_jacobian_policy_regular_stepper<
  T, ::pressio::ode::implicitmethods::Euler, args...>;

template<typename T, typename ... args>
using admissible_implicit_bdf2_jacobian_policy_regular_stepper =
  admissible_implicit_jacobian_policy_regular_stepper<
  T, ::pressio::ode::implicitmethods::BDF2, args...>;

}}} // namespace pressio::ode::meta
#endif
