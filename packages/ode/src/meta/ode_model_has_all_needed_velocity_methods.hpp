
#ifndef ODE_MODEL_HAS_ALL_NEEDED_VELOCITY_METHODS_HPP_
#define ODE_MODEL_HAS_ALL_NEEDED_VELOCITY_METHODS_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../mpl/src/detection_idiom.hpp"

namespace rompp{ namespace ode{ namespace meta {


template <typename T, typename a_t, typename b_t, typename = void>
struct has_velocity_method_callable_with_two_args : std::false_type{};

template <typename T, typename a_t, typename b_t>
struct has_velocity_method_callable_with_two_args<
  T, a_t, b_t,
  ::rompp::mpl::enable_if_t<
    !std::is_void<
      decltype(
	       std::declval<T>().velocity(
            std::declval<a_t const&>(),
					  std::declval<b_t &>())
	   )
      >::value
    >
  > : std::true_type{};



template <typename T, typename a_t,
	  typename b_t, typename c_t, typename = void>
struct has_velocity_method_callable_with_three_args : std::false_type{};

template <typename T, typename a_t, typename b_t, typename c_t>
struct has_velocity_method_callable_with_three_args<
  T, a_t, b_t, c_t,
  ::rompp::mpl::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T>().velocity(
					  std::declval<a_t const&>(),
					  std::declval<b_t &>(),
					  std::declval<c_t &>()
					  )
	   )
      >::value
    >
  > : std::true_type{};

//---------------------------------------------------------------


template<typename model_type,
	 typename state_type,
	 typename velocity_type,
	 typename scalar_type,
	 typename enable = void>
struct model_has_needed_velocity_methods : std::false_type{};

template<typename model_type,
	 typename state_type,
	 typename velocity_type,
	 typename scalar_type>
struct model_has_needed_velocity_methods<
  model_type, state_type, velocity_type, scalar_type,
  typename std::enable_if<
   // has method with 2 arguments,
    has_velocity_method_callable_with_two_args<
     model_type, state_type, scalar_type
     >::value and
    // containers::meta::is_detected<
    // has_velocity_method_callable_with_two_args,
    //  model_type, state_type, scalar_type
    //  >::value and
   // has velocity method with 3 arguments
    /*containers::meta::is_detected<*/
     has_velocity_method_callable_with_three_args<
     model_type, state_type, velocity_type, scalar_type
     >::value and
   // method with 2 arguments returns a velocity_type
   std::is_same<
      velocity_type,
      decltype(
	       std::declval<model_type>().velocity
	       (std::declval<state_type const&>(),
		std::declval<scalar_type &>())
	   )
      >::value
  >::type
  > : std::true_type{};


}}} // namespace rompp::ode::meta
#endif
