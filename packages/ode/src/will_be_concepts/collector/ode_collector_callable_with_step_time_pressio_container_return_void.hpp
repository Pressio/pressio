
#ifndef ode_collector_callable_with_step_time_pressio_container_HPP_
#define ode_collector_callable_with_step_time_pressio_container_HPP_

namespace pressio{ namespace ode{ namespace meta {

template<
  typename T,
  typename time_type,
  typename state_type,
  typename enable = void
  >
struct collector_callable_with_step_time_pressio_container_return_void
  : std::false_type{};


template<
  typename T,
  typename time_type,
  typename state_type
  >
struct collector_callable_with_step_time_pressio_container_return_void<
  T, time_type, state_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_wrapper<state_type>::value and
    std::is_void<
      decltype
      (
       std::declval<T>()
       (
       std::declval<::pressio::ode::types::step_t>(),
       std::declval<time_type>(),
       std::declval<const state_type &>()
	     )
       )
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::ode::meta
#endif
