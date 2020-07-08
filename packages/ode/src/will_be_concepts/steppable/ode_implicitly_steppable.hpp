
#ifndef ode_implicitly_steppable_HPP_
#define ode_implicitly_steppable_HPP_

namespace pressio{ namespace ode{ namespace concepts {

template <
  typename T, typename state_type, typename time_type, typename solver_type,
  typename enable = void
  >
struct implicitly_steppable
  : std::false_type{};

template <
  typename T, typename state_type, typename time_type, typename solver_type
  >
struct implicitly_steppable<
  T, state_type, time_type, solver_type,
  mpl::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T>().doStep
              (
              std::declval<state_type &>(),
              std::declval<time_type const &>(),
              std::declval<time_type const &>(),
						  std::declval<::pressio::ode::types::step_t const &>(),
              std::declval<solver_type &>()
						  )
	       )
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::ode::concepts
#endif
