
#ifndef ODE_MODEL_HAS_ALL_NEEDED_VELOCITY_METHODS_HPP_
#define ODE_MODEL_HAS_ALL_NEEDED_VELOCITY_METHODS_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../mpl/src/detection_idiom.hpp"

namespace pressio{ namespace ode{ namespace meta {


template <typename T, typename state_t, typename sc_t, typename = void>
struct has_velocity_method_callable_with_two_args : std::false_type{};

template <typename T, typename state_t, typename sc_t>
struct has_velocity_method_callable_with_two_args<
  T, state_t, sc_t,
  ::pressio::mpl::enable_if_t<
    !std::is_void<
      decltype(
	       std::declval<T>().velocity(
					  std::declval<state_t const&>(),
					  std::declval<sc_t>()
					  )
	       )
      >::value
    >
  > : std::true_type{};



template <typename T,
	  typename state_t,
	  typename sc_t,
	  typename velo_t,
	  typename = void>
struct has_velocity_method_callable_with_three_args : std::false_type{};

template <typename T,
	  typename state_t,
	  typename sc_t,
	  typename velo_t>
struct has_velocity_method_callable_with_three_args<
  T, state_t, sc_t, velo_t,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T>().velocity(
					  std::declval<state_t const&>(),
					  std::declval<sc_t>(),
					  std::declval<velo_t &>()
					  )
	   )
      >::value
    >
  > : std::true_type{};

//---------------------------------------------------------------


template<
  typename model_type,
  typename state_type,
  typename velocity_type,
  typename scalar_type,
  typename enable = void
  >
struct model_has_needed_velocity_methods : std::false_type{};

template<
  typename model_type,
  typename state_type,
  typename velocity_type,
  typename scalar_type
  >
struct model_has_needed_velocity_methods<
  model_type, state_type, velocity_type, scalar_type,
  mpl::enable_if_t<
    // has method with 2 arguments,
    has_velocity_method_callable_with_two_args<
      model_type, state_type, scalar_type
      >::value and
    // has velocity method with 3 arguments
    has_velocity_method_callable_with_three_args<
      model_type, state_type, scalar_type, velocity_type
      >::value and
    // method with 2 arguments returns a velocity_type
    mpl::is_same<
      velocity_type,
      decltype(
	       std::declval<model_type>().velocity
	       (
		std::declval<state_type const&>(),
		std::declval<scalar_type>())
	       )
      >::value
    >
  > : std::true_type{};


}}} // namespace pressio::ode::meta
#endif
