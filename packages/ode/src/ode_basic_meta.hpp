
#ifndef ODE_META_BASIC_META_HPP_
#define ODE_META_BASIC_META_HPP_

#include "ode_ConfigDefs.hpp"
#include "../../core/src/vector/core_vector_meta.hpp"
#include "../../core/src/multi_vector/core_multi_vector_meta.hpp"
#include "../../core/src/matrix/core_matrix_meta.hpp"

namespace rompp{ namespace ode{ namespace meta {

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
template<typename model_type,
	 typename enable = void>
struct is_legitimate_model_for_explicit_ode : std::false_type{};

template<typename model_type>
struct is_legitimate_model_for_explicit_ode<model_type,
  core::meta::enable_if_t<
    !core::meta::is_core_vector_wrapper<model_type>::value
    >>
  : std::true_type{};
      
//-------------------------------------------------------
template<typename model_type,
	 typename enable = void>
struct is_legitimate_model_for_implicit_ode : std::false_type{};

template<typename model_type>
struct is_legitimate_model_for_implicit_ode<model_type,
  core::meta::enable_if_t<
    !core::meta::is_core_vector_wrapper<model_type>::value
    >>
  : std::true_type{};
      
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
      
 
} // namespace meta
} // namespace ode
} //end namespace rompp
#endif
