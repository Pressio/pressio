
#ifndef ODE_META_META_HPP_
#define ODE_META_META_HPP_

#include "../../../core/src/meta/core_meta_basic.hpp"

namespace ode{
namespace meta {

// template<typename time_type, typename enable = void>
// struct isLegitimateTimeType : std::false_type{};

// template<typename time_type>
// struct isLegitimateTimeType<
//   time_type,
//   typename std::enable_if<
//     std::is_floating_point<time_type>::value
//     >::type
//   > : std::true_type{};

//---------------------------------------------------------------

// when providing a collector functor for the integrator,
// this has to be a functor. With following syntax:
//
//  void operator()(step, time, state)
//
    
template<typename functor,
	 typename int_type,
	 typename time_type,
	 typename state_type,
	 typename enable = void>
struct isLegitimateCollector : std::false_type{};
  
template<typename functor,
	 typename int_type,
	 typename time_type,
	 typename state_type>
struct isLegitimateCollector<
  functor, int_type, time_type, state_type,
  core::meta::void_t<
    decltype(std::declval<functor>()(
      std::declval<int_type>(),
      std::declval<time_type>(),
      std::declval<state_type>()
	      )
	    )
    >
  > : std::true_type{};
  
 
} // namespace meta
} // namespace core
#endif
