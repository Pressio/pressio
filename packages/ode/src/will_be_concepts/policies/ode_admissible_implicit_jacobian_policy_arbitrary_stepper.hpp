
#ifndef ode_admissible_jacobian_policy_for_implicit_arbitrary_stepper_HPP_
#define ode_admissible_jacobian_policy_for_implicit_arbitrary_stepper_HPP_

namespace pressio{ namespace ode{ namespace meta {

template<
  typename T,
  typename state_t,
  typename jacobian_t,
  typename system_t,
  typename enable = void
  >
struct jacobian_policy_callable_with_four_args : std::false_type{};

template<
  typename T,
  typename state_t,
  typename jacobian_t,
  typename system_t
  >
struct jacobian_policy_callable_with_four_args<
  T, state_t, jacobian_t, system_t,
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
    >
  > : std::true_type{};
/////////////////////////

template<
  typename T,
  std::size_t numPrevStates,
  typename state_t,
  typename jacobian_t,
  typename system_t,
  typename scalar_t,
  typename enable = void
  >
struct jacobian_policy_callable_with_seven_args : std::false_type{};

template<
  typename T,
  std::size_t numPrevStates,
  typename state_t,
  typename jacobian_t,
  typename system_t,
  typename scalar_t
  >
struct jacobian_policy_callable_with_seven_args<
  T, numPrevStates, state_t, jacobian_t, system_t, scalar_t,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const>().template operator()
       (
	std::declval<state_t const &>(),
	std::declval<::pressio::ode::AuxStatesContainer<false, state_t, numPrevStates> const &>(),
	std::declval<system_t const &>(),
	std::declval<scalar_t const & >(),
	std::declval<scalar_t const &>(),
	std::declval<::pressio::ode::types::step_t const &>(),
	std::declval<jacobian_t &>()
	)
       )
      >::value
    >
  > : std::true_type{};
/////////////////////////


template<
  typename T,
  std::size_t numPrevStates,
  typename state_t,
  typename jacobian_t,
  typename system_t,
  typename scalar_t,
  typename enable = void
  >
struct admissible_jacobian_policy_for_implicit_arbitrary_stepper
{
  static constexpr auto c1 = ::pressio::ode::meta::legitimate_implicit_state_type<state_t>::value;
  static constexpr auto c2 = ::pressio::ode::meta::legitimate_jacobian_type<jacobian_t>::value;

  static constexpr auto c5 = jacobian_policy_callable_with_four_args<
    T, state_t, jacobian_t, system_t>::value;

  static constexpr auto c6 = jacobian_policy_callable_with_seven_args<
    T, numPrevStates, state_t, jacobian_t, system_t, scalar_t>::value;

  using value_type = bool;
  static constexpr value_type value = c1 && c2 && c5 && c6;
  using type = std::integral_constant<value_type, value>;
};

}}} // namespace pressio::ode::meta
#endif
