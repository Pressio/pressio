
#ifndef ode_explicitly_steppable_HPP_
#define ode_explicitly_steppable_HPP_

namespace pressio{ namespace ode{ namespace concepts {

template <
  typename T, typename state_type, typename time_type, 
  typename enable = void
  >
struct explicitly_steppable
  : std::false_type{};

template <
  typename T, typename state_type, typename time_type
  >
struct explicitly_steppable<
  T, state_type, time_type, 
  mpl::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T>().doStep(
              std::declval<state_type &>(),
              std::declval<time_type const &>(),
              std::declval<time_type const &>(),
						  std::declval<::pressio::ode::types::step_t const &>()
						  )
	       )
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::ode::concepts
#endif
