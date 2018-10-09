
#ifndef ODE_MODEL_HAS_ALL_NEEDED_RESIDUAL_METHODS_HPP_
#define ODE_MODEL_HAS_ALL_NEEDED_RESIDUAL_METHODS_HPP_

#include "ode_ConfigDefs.hpp"

namespace rompp{ namespace ode{ namespace meta {
      
template <typename T, typename a_t, typename b_t>
using has_residual_method_callable_with_two_args =
  decltype(std::declval<T>().residual(std::declval<a_t const&>(),
				      std::declval<b_t &>())
	   );

template <typename T,
	  typename a_t, typename b_t, typename c_t>
using has_residual_method_callable_with_three_args =
  decltype(std::declval<T>().residual(std::declval<a_t const&>(),
				      std::declval<b_t &>(),
				      std::declval<c_t>())
	   );
      
//-------------------------------------------------------

template<typename model_type,
	 typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename enable = void>
struct model_has_needed_residual_methods : std::false_type{};

template<typename model_type,
	 typename state_type,
	 typename residual_type,
	 typename scalar_type>
struct model_has_needed_residual_methods<
  model_type, state_type, residual_type, scalar_type,
  typename std::enable_if<
   // has residual method with 3 arguments
   core::meta::is_detected<
     has_residual_method_callable_with_three_args,
     model_type, state_type, residual_type, scalar_type
     >::value and 
   // has residual method with 2 arguments,
   core::meta::is_detected<
    has_residual_method_callable_with_two_args,
     model_type, state_type, scalar_type
     >::value and
   // residual method with 2 arguments returns a residual_type
   std::is_same<
    core::meta::detected_t<
      has_residual_method_callable_with_two_args,
      model_type, state_type, scalar_type>,
     residual_type
    >::value    
  >::type
  > : std::true_type{};

//------------------------------------------------------
       
}}} // namespace rompp::ode::meta
#endif
