
#ifndef ODE_COLLECTOR_ACCEPTS_NATIVE_CONTAINER_HPP_
#define ODE_COLLECTOR_ACCEPTS_NATIVE_CONTAINER_HPP_

#include "../ode_fwd.hpp"

namespace pressio{ namespace ode{ namespace meta {

template<
  typename collector_type,
  typename int_type,
  typename time_type,
  typename state_type,
  typename enable = void
  >
struct collector_accepts_native_container
  : std::false_type{};


/* true if state_type is a wrapper and collector accepts the wrapped type */
template<
  typename collector_type,
  typename int_type,
  typename time_type,
  typename state_type
  >
struct collector_accepts_native_container<
  collector_type, int_type, time_type, state_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_wrapper<state_type>::value and
    std::is_void<
      decltype
      (
       std::declval<collector_type>()
       (
	std::declval<int_type>(),
	std::declval<time_type>(),
	std::declval<
	const typename ::pressio::containers::details::traits<state_type>::wrapped_t &
	>()
	)
       )
      >::value
    >
  > : std::true_type{};


/* true if state_type is NOT wrapper and collector accepts the type itself */
template<
  typename collector_type,
  typename int_type,
  typename time_type,
  typename state_type
  >
struct collector_accepts_native_container<
  collector_type, int_type, time_type, state_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_wrapper<state_type>::value == false and
    std::is_void<
      decltype
      (
       std::declval<collector_type>()
       (
	std::declval<int_type>(),
	std::declval<time_type>(),
	std::declval<const state_type &>()
	)
       )
      >::value
    >
  > : std::true_type{};


}}} // namespace pressio::ode::meta
#endif
