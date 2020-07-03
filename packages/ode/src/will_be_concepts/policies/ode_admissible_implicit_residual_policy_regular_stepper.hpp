
#ifndef ODE_ode_admissible_implicit_residual_policy_regular_stepper_HPP_
#define ODE_ode_admissible_implicit_residual_policy_regular_stepper_HPP_

namespace pressio{ namespace ode{ namespace meta {

template<
  typename T,
  typename tag,
  std::size_t numPrevStates,
  typename state_t,
  typename residual_t,
  typename system_t,
  typename scalar_t,
  typename enable = void
  >
struct admissible_implicit_residual_policy_regular_stepper 
  : std::false_type{};


template<
  typename T,
  typename tag,
  std::size_t numPrevStates,
  typename state_t,
  typename residual_t,
  typename system_t,
  typename scalar_t
  >
struct admissible_implicit_residual_policy_regular_stepper<
  T, tag, numPrevStates, state_t, residual_t, system_t, scalar_t,
  ::pressio::mpl::enable_if_t<
    // is callable with state and system
    std::is_same<
      residual_t,
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
    // is callable with state, system, etc
    std::is_void<
      decltype
      (
       std::declval<T const>().template operator()
       <tag>(
	     std::declval<state_t const &>(),
	     std::declval<::pressio::ode::AuxStatesContainer<false, state_t, numPrevStates> const &>(),
	     std::declval<system_t const &>(),
	     std::declval<scalar_t const &>(),
	     std::declval<scalar_t const &>(),
	     std::declval<::pressio::ode::types::step_t>(),
	     std::declval<residual_t &>(),
	     ::pressio::solvers::Norm::Undefined,
	     std::declval<scalar_t &>()
	     )
       )
      >::value
    >
  > : std::true_type{};
//------------------------------------------------------------------


template<typename T, typename ... args>
using admissible_implicit_euler_residual_policy_regular_stepper =
  admissible_implicit_residual_policy_regular_stepper<
  T, ::pressio::ode::implicitmethods::Euler, 1, args...>;

template<typename T, typename ... args>
using admissible_implicit_bdf2_residual_policy_regular_stepper =
  admissible_implicit_residual_policy_regular_stepper<
  T, ::pressio::ode::implicitmethods::BDF2, 2, args...>;

}}} // namespace pressio::ode::meta
#endif
