
#ifndef ODE_IS_LEGITIMATE_COLLECTOR_HPP_
#define ODE_IS_LEGITIMATE_COLLECTOR_HPP_

#include "ode_collector_accepts_native_container.hpp"
#include "ode_collector_accepts_pressio_container.hpp"

namespace pressio{ namespace ode{ namespace meta {

template<
  typename collector_type,
  typename int_type,
  typename time_type,
  typename state_type
  >
struct is_legitimate_collector{

  static constexpr auto collector_accepting_native_container =
    ::pressio::ode::meta::collector_accepts_native_container<collector_type, int_type, time_type, state_type>::value;

  static constexpr auto collector_accepting_pressio_container =
    ::pressio::ode::meta::collector_accepts_pressio_container<collector_type, int_type, time_type, state_type>::value;

  // force user to only use one or the other for now
  static constexpr auto both_are_true = (collector_accepting_native_container and collector_accepting_pressio_container);
  static_assert( both_are_true == false,
		 "Currently, the collector/observer passed to ode must either accept a native container or a pressio container wrapper. You cannot have two methods to cover both cases. Pick one. ");

  // value is true if either one is true
  static constexpr auto value = collector_accepting_native_container or collector_accepting_pressio_container;
};

}}} // namespace pressio::ode::meta
#endif
