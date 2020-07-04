

#ifndef ODE_ode_admissible_residual_policy_for_implicit_arbitrary_stepper_HPP_
#define ODE_ode_admissible_residual_policy_for_implicit_arbitrary_stepper_HPP_

namespace pressio{ namespace ode{ namespace meta {

template<
  typename T,
  typename residual_t,
  typename system_t,
  typename enable = void
  >
struct residual_policy_create_method : std::false_type{};

template<
  typename T,
  typename residual_t,
  typename system_t
  >
struct residual_policy_create_method<
  T, residual_t, system_t,
  ::pressio::mpl::enable_if_t<
    std::is_same<
      residual_t,
      decltype
      (
       std::declval<T const>().create(std::declval<system_t const &>())
       )
      >::value
    >
  > : std::true_type{};
/////////////////////////

template<
  typename T,
  std::size_t numPrevStates,
  typename state_t,
  typename residual_t,
  typename system_t,
  typename scalar_t,
  typename enable = void
  >
struct residual_policy_compute_method : std::false_type{};

template<
  typename T,
  std::size_t numPrevStates,
  typename state_t,
  typename residual_t,
  typename system_t,
  typename scalar_t
  >
struct residual_policy_compute_method<
  T, numPrevStates, state_t, residual_t, system_t, scalar_t,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype
      (
       std::declval<T const>().template compute
       (
	std::declval<state_t const &>(),
	std::declval<::pressio::ode::AuxStatesContainer<false, state_t, numPrevStates> const &>(),
	std::declval<system_t const &>(),
	std::declval<scalar_t const &>(),
	std::declval<scalar_t const &>(),
	std::declval<::pressio::ode::types::step_t const &>(),
	std::declval<residual_t &>(),
	::pressio::Norm::Undefined,
	std::declval<scalar_t &>()
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
  typename residual_t,
  typename system_t,
  typename scalar_t,
  typename enable = void
  >
struct admissible_residual_policy_for_implicit_arbitrary_stepper
{
  static constexpr bool c1 = ::pressio::ode::meta::legitimate_implicit_state_type<state_t>::value;
  static constexpr bool c2 = ::pressio::ode::meta::legitimate_implicit_residual_type<residual_t>::value;

  static constexpr bool c5 = residual_policy_create_method<T, residual_t, system_t>::value;

  static constexpr bool c6 = residual_policy_compute_method<
    T, numPrevStates, state_t, residual_t, system_t, scalar_t>::value;

  using value_type = bool;
  static constexpr value_type value = c1 && c2 && c5 && c6;
  using type = std::integral_constant<value_type, value>;
};

}}} // namespace pressio::ode::meta
#endif
