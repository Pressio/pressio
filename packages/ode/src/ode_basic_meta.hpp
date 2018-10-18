
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
      
//---------------------------------------------------------
template<typename state_type, typename enable = void>
struct is_legitimate_explicit_state_type : std::false_type{};

template<typename state_type>
struct is_legitimate_explicit_state_type<state_type,
 typename std::enable_if<
   core::meta::is_core_vector_wrapper<state_type>::value
   >::type
  > : std::true_type{};

      
//---------------------------------------------------------
template<typename residual_type, typename enable = void>
struct is_legitimate_explicit_residual_type : std::false_type{};

template<typename residual_type>
struct is_legitimate_explicit_residual_type<residual_type,
 typename std::enable_if<
   core::meta::is_core_vector_wrapper<residual_type>::value
   >::type
  > : std::true_type{};
      

//---------------------------------------------------------
template<typename state_type, typename enable = void>
struct is_legitimate_implicit_state_type : std::false_type{};

template<typename state_type>
struct is_legitimate_implicit_state_type<state_type,
        typename std::enable_if<
	  core::meta::is_core_vector_wrapper<state_type>::value
	  >::type > : std::true_type{};

      
//---------------------------------------------------------
template<typename residual_type, typename enable = void>
struct is_legitimate_implicit_residual_type : std::false_type{};

template<typename residual_type>
struct is_legitimate_implicit_residual_type<residual_type,
 typename std::enable_if<
   core::meta::is_core_vector_wrapper<residual_type>::value
   >::type
  > : std::true_type{};

      
//---------------------------------------------------------
template<typename jacobian_type, typename enable = void>
struct is_legitimate_jacobian_type : std::false_type{};

template<typename jacobian_type>
struct is_legitimate_jacobian_type<jacobian_type,
       typename std::enable_if<
	 core::meta::is_core_matrix_wrapper<jacobian_type>::value or
	 core::meta::is_core_multi_vector_wrapper<jacobian_type>::value
	 >::type
       > : std::true_type{};
 
}}} // namespace rompp::ode::meta
#endif
