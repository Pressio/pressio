
#ifndef ODE_IS_LEGITIMATE_COLLECTOR_HPP_
#define ODE_IS_LEGITIMATE_COLLECTOR_HPP_

#include "ode_basic_meta.hpp"
#include "ode_model_has_all_needed_residual_methods.hpp"

namespace rompp{ namespace ode{ namespace meta {

/*-------------------------------------------------------
 * when providing a collector functor for the integrator,
 * this has to be a functor. With following syntax:
 *  void operator()(step, time, state)
 */
template<typename collector_type,
	 typename int_type,
	 typename time_type,
	 typename state_type,
	 typename enable = void>
struct is_legitimate_collector : std::false_type{};

template<typename collector_type,
	 typename int_type,
	 typename time_type,
	 typename state_type>
struct is_legitimate_collector<
  collector_type, int_type, time_type, state_type,
  ::rompp::mpl::void_t<
    decltype(std::declval<collector_type>()(std::declval<int_type>(),
					    std::declval<time_type>(),
					    std::declval<const state_type &>()
					    )
	     )
    >
  > : std::true_type{};


}}} // namespace rompp::ode::meta
#endif
