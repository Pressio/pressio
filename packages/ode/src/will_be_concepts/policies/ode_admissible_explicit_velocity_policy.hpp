
#ifndef ODE_ode_admissible_explicit_velocity_policy_HPP_
#define ODE_ode_admissible_explicit_velocity_policy_HPP_

namespace pressio{ namespace ode{ namespace meta {

template<
  typename T,
  typename scalar_t, typename state_t, typename velocity_t, typename system_t,
  typename enable = void
  >
struct admissible_explicit_velocity_policy
  : std::false_type{};

template<
  typename T,
  typename scalar_t, typename state_t, typename velocity_t, typename system_t
  >
struct admissible_explicit_velocity_policy<
  T, scalar_t, state_t, velocity_t, system_t,
  mpl::enable_if_t<
    std::is_same<
      decltype
      (
       std::declval<T const>().create(std::declval<system_t const &>())
       ),
      velocity_t
      >::value
    and
    std::is_void<
      decltype
      (
       std::declval<T const>().template compute
       (
      	std::declval<state_t const &>(),
      	std::declval<velocity_t &>(),
      	std::declval<system_t const &>(),
      	std::declval<scalar_t const &>()
      	)
       )
      >::value
    >
  > : std::true_type{};


}}}//end namespace pressio::containers::meta
#endif
