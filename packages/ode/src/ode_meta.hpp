
#ifndef ODE_META_HPP_
#define ODE_META_HPP_

#include <type_traits>
#include "meta/core_meta_basic.hpp"

namespace ode{
namespace meta {

template<typename state_type, typename enable = void>
struct isLegitimateStateType : std::false_type{};

template<typename state_type>
struct isLegitimateStateType<
  state_type,
  typename std::enable_if<!std::is_integral<state_type>::value &&
			  !std::is_floating_point<state_type>::value
			  >::type
  > : std::true_type{};
//---------------------------------------------------------------

template<typename residual_type, typename enable = void>
struct isLegitimateResidualType : std::false_type{};

template<typename residual_type>
struct isLegitimateResidualType<
  residual_type,
  typename std::enable_if<!std::is_integral<residual_type>::value &&
			  !std::is_floating_point<residual_type>::value
			  >::type
  > : std::true_type{};
//---------------------------------------------------------------


template<typename time_type, typename enable = void>
struct isLegitimateTimeType : std::false_type{};

template<typename time_type>
struct isLegitimateTimeType<
  time_type,
  typename std::enable_if<std::is_floating_point<time_type>::value>::type
  > : std::true_type{};
//---------------------------------------------------------------
  

} // namespace meta
} // namespace core
#endif
