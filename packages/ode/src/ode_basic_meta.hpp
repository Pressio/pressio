
#ifndef ODE_META_BASIC_META_HPP_
#define ODE_META_BASIC_META_HPP_

#include "ode_ConfigDefs.hpp"
#include "../../core/src/meta/core_meta_detection_idiom.hpp"
#include "../../core/src/vector/core_vector_meta.hpp"
#include "../../core/src/multi_vector/core_multi_vector_meta.hpp"
#include "../../core/src/matrix/core_matrix_meta.hpp"

namespace rompp{ namespace ode{ namespace meta {

template <typename T>
using has_state_typedef = typename T::state_type;

template <typename T>
using has_residual_typedef = typename T::residual_type;

template <typename T>
using has_jacobian_typedef = typename T::jacobian_type;

template <typename T>
using has_scalar_typedef = typename T::scalar_type;
      

//---------------------------------------------------------------
template<typename state_type, typename enable = void>
struct is_legitimate_explicit_state_type : std::false_type{};

template<typename state_type>
struct is_legitimate_explicit_state_type<state_type,
 typename std::enable_if<
   core::meta::is_core_vector_wrapper<state_type>::value
   >::type
  > : std::true_type{};
  
//---------------------------------------------------------------
// residual satisfies same constraints as state type (for now)
template<typename residual_type, typename enable = void>
struct is_legitimate_explicit_residual_type
  : is_legitimate_explicit_state_type<residual_type>{};


//--------------------------------------------------------------
template<typename state_type, typename enable = void>
struct is_legitimate_implicit_state_type : std::false_type{};

template<typename state_type>
struct is_legitimate_implicit_state_type<state_type,
        typename std::enable_if<
	  core::meta::is_core_vector_wrapper<state_type>::value
	  >::type > : std::true_type{};

//---------------------------------------------------------------
// residual satisfies same constraints as state type (for now)
template<typename residual_type, typename enable = void>
struct is_legitimate_implicit_residual_type
  : is_legitimate_implicit_state_type<residual_type>{};

//---------------------------------------------------------------
template<typename jacobian_type, typename enable = void>
struct is_legitimate_jacobian_type : std::false_type{};

template<typename jacobian_type>
struct is_legitimate_jacobian_type<jacobian_type,
       typename std::enable_if<
	 core::meta::is_core_matrix_wrapper<jacobian_type>::value
	 >::type
       > : std::true_type{};


      
//-------------------------------------------------------
// when providing a collector functor for the integrator,
// this has to be a functor. With following syntax:
//
//  void operator()(step, time, state)  
template<typename functor,
	 typename int_type,
	 typename time_type,
	 typename state_type,
	 typename enable = void>
struct is_legitimate_collector : std::false_type{};
  
template<typename functor,
	 typename int_type,
	 typename time_type,
	 typename state_type>
struct is_legitimate_collector<
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

//-------------------------------------------------------
 
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
	 typename enable = void>
struct model_has_needed_residual_methods : std::false_type{};

template<typename model_type>
struct model_has_needed_residual_methods<
  model_type,
  typename std::enable_if<
   // has residual method with 3 arguments
   core::meta::is_detected<
     has_residual_method_callable_with_three_args,
     model_type,
     typename model_type::state_type,
     typename model_type::residual_type,
     typename model_type::scalar_type>::value and 
   // has residual method with 2 arguments,
   core::meta::is_detected<
    has_residual_method_callable_with_two_args,
     model_type,
     typename model_type::state_type,
     typename model_type::scalar_type
   >::value and
   // residual method with 2 arguments returns a residual_type
   std::is_same<
    core::meta::detected_t<
      has_residual_method_callable_with_two_args,
      model_type, typename model_type::state_type,
      typename model_type::scalar_type>,
     typename model_type::residual_type
    >::value    
  >::type
  > : std::true_type{};

//------------------------------------------------------

      
template<typename model_type,
	 typename enable = void>
struct is_legitimate_model_for_explicit_ode : std::false_type{};

template<typename model_type>
struct is_legitimate_model_for_explicit_ode<
  model_type,
  typename std::enable_if<
   // has to have scalar typedef
   core::meta::is_detected<has_scalar_typedef, model_type>::value and 
   // has to have state typedef
   core::meta::is_detected<has_state_typedef, model_type>::value and 
   // has to have residual typedef
   core::meta::is_detected<has_residual_typedef, model_type>::value and 
   // has residual methods
   model_has_needed_residual_methods<model_type>::value
  >::type
  > : std::true_type{};

//-------------------------------------------------------



template <typename T, typename a_t, typename b_t>
using has_jacobian_method_callable_with_two_args =
  decltype(std::declval<T>().jacobian(std::declval<a_t const&>(),
				      std::declval<b_t &>())
	   );

template <typename T,
	  typename a_t, typename b_t, typename c_t>
using has_jacobian_method_callable_with_three_args =
  decltype(std::declval<T>().jacobian(std::declval<a_t const&>(),
				      std::declval<b_t &>(),
				      std::declval<c_t>())
	   );
      
//-------------------------------------------------------

      
template<typename model_type,
	 typename enable = void>
struct model_has_needed_jacobian_methods : std::false_type{};

template<typename model_type>
struct model_has_needed_jacobian_methods<
  model_type,
  typename std::enable_if<
   // has jacobian method with 3 arguments
   core::meta::is_detected<
     has_jacobian_method_callable_with_three_args,
     model_type,
     typename model_type::state_type,
     typename model_type::jacobian_type,
     typename model_type::scalar_type>::value and 
   // has jacobian method with 2 arguments,
   core::meta::is_detected<
    has_jacobian_method_callable_with_two_args,
     model_type,
     typename model_type::state_type,
     typename model_type::scalar_type
   >::value and
   // jacobian method with 2 arguments returns a residual_type
   std::is_same<
    core::meta::detected_t<
      has_jacobian_method_callable_with_two_args,
      model_type,
      typename model_type::state_type,
      typename model_type::scalar_type>,
      typename model_type::jacobian_type
    >::value    
  >::type
  > : std::true_type{};

//------------------------------------------------------
      
template<typename model_type,
	 typename enable = void>
struct is_legitimate_model_for_implicit_ode : std::false_type{};

template<typename model_type>
struct is_legitimate_model_for_implicit_ode<
  model_type,
  typename std::enable_if<
   // has to have scalar typedef
   core::meta::is_detected<has_scalar_typedef, model_type>::value and 
   // has to have state typedef
   core::meta::is_detected<has_state_typedef, model_type>::value and 
   // has to have residual typedef
   core::meta::is_detected<has_residual_typedef, model_type>::value and 
   // has to have jacobian typedef
   core::meta::is_detected<has_jacobian_typedef, model_type>::value and 
   // has residual methods
   model_has_needed_residual_methods<model_type>::value
  >::type
  > : std::true_type{};


//------------------------------------------------------

// for now, just leave Backward Euler as the legitimate stepper
template<typename T,
	 typename enable = void>
struct is_legitimate_auxiliary_stepper : std::false_type{};

template<typename T>
struct is_legitimate_auxiliary_stepper<T,
  typename std::enable_if<
    ::rompp::ode::details::traits<typename T::base_t>::enum_id ==
    ::rompp::ode::ImplicitSteppersEnum::Euler
    >::type
  > : std::true_type{};
      
      
 
}}} // namespace rompp::ode::meta
#endif
