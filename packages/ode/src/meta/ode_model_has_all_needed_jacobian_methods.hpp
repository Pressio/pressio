
#ifndef ODE_MODEL_HAS_ALL_NEEDED_JACOBIAN_METHODS_HPP_
#define ODE_MODEL_HAS_ALL_NEEDED_JACOBIAN_METHODS_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../mpl/src/detection_idiom.hpp"

namespace pressio{ namespace ode{ namespace meta {

template <typename T, typename state_t, typename sc_t, typename = void>
struct has_jacobian_method_callable_with_two_args : std::false_type{};

template <typename T, typename state_t, typename sc_t>
struct has_jacobian_method_callable_with_two_args<
  T, state_t, sc_t,
  ::pressio::mpl::void_t<
  decltype(
	   std::declval<T>().jacobian
	   (std::declval<state_t const&>(),
	    std::declval<sc_t>())
	   )
    >
  > : std::true_type{};

//-------------------------------------------------------

template <typename T,
	  typename state_t,
	  typename sc_t,
	  typename jac_t,
	  typename = void>
struct has_jacobian_method_callable_with_three_args : std::false_type{};

template <typename T,
	  typename state_t,
	  typename sc_t,
	  typename jac_t>
struct has_jacobian_method_callable_with_three_args<
  T, state_t, sc_t, jac_t,
  ::pressio::mpl::void_t<
  decltype(
	   std::declval<T>().jacobian
	   (std::declval<state_t const&>(),
	    std::declval<sc_t>(),
	    std::declval<jac_t &>())
	   )
    >
  >: std::true_type{};

//-------------------------------------------------------


template<typename model_type,
	 typename state_type,
	 typename jacobian_type,
	 typename scalar_type,
	 typename enable = void>
struct model_has_needed_jacobian_methods : std::false_type{};

template<typename model_type,
	 typename state_type,
	 typename jacobian_type,
	 typename scalar_type>
struct model_has_needed_jacobian_methods<
  model_type, state_type, jacobian_type, scalar_type,
  typename std::enable_if<
    // has jacobian method with 2 arguments,
    has_jacobian_method_callable_with_two_args<
      model_type, state_type, scalar_type
      >::value and
    // has jacobian method with 3 arguments
    has_jacobian_method_callable_with_three_args<
      model_type, state_type, scalar_type, jacobian_type
      >::value and
    // jacobian method with 2 arguments returns a residual_type
    std::is_same<
      jacobian_type,
      decltype(std::declval<model_type>().jacobian
	       (
		std::declval<state_type const&>(),
		std::declval<scalar_type>()
		)
	       )
      >::value
  >::type
  > : std::true_type{};

//------------------------------------------------------

}}} // namespace pressio::ode::meta
#endif
