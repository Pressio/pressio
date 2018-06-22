
#ifndef ODE_META_HPP_
#define ODE_META_HPP_

#include <type_traits>
#include "meta/core_meta_basic.hpp"
#include "vector/core_vector_traits.hpp"
#include "matrix/core_matrix_traits.hpp"


namespace ode{
namespace meta {


template<typename time_type, typename enable = void>
struct isLegitimateTimeType : std::false_type{};

template<typename time_type>
struct isLegitimateTimeType<
  time_type,
  typename std::enable_if<std::is_floating_point<time_type>::value>::type
  > : std::true_type{};


//********************************************
//         FOR EXPLICIT METHODS
//********************************************
  
// For explicit methods, things are easier so the admissible types
// for states is kind of anything as long as has [] operator.
// Because we just wrap it inside, without needing to do much else.
template<typename state_type, typename enable = void>
struct isLegitimateExplicitStateType : std::false_type{};

template<typename state_type>
struct isLegitimateExplicitStateType<state_type,
 typename std::enable_if<core::meta::is_coreVectorWrapper<state_type>::value>::type
  > : std::true_type{};

template<typename state_type>
struct isLegitimateExplicitStateType<state_type,
  typename std::enable_if<core::meta::is_vectorStdLib<state_type>::value>::type
  > : std::true_type{};

template<typename state_type>
struct isLegitimateExplicitStateType<state_type,
  typename std::enable_if<core::meta::is_vectorEigen<state_type>::value>::type
  > : std::true_type{};

template<typename state_type>
struct isLegitimateExplicitStateType<state_type,
  typename std::enable_if<core::meta::is_vectorEpetra<state_type>::value>::type
  > : std::true_type{};

//---------------------------------------------------------------
// residual has to be same type as state for now

template<typename residual_type, typename enable = void>
struct isLegitimateExplicitResidualType
  : isLegitimateExplicitStateType<residual_type>{};



//********************************************
//         FOR IMPLICIT METHODS
//********************************************
// For implicit methods, things are more complicated. We need to solve
// linear or non linear systems. So allow only if the state is wrapped
// with out core::vector class. 

template<typename state_type, typename enable = void>
struct isLegitimateImplicitStateType : std::false_type{};

template<typename state_type>
struct isLegitimateImplicitStateType<state_type,
        typename std::enable_if<
	  core::meta::is_coreVectorWrapper<state_type>::value
	  >::type > : std::true_type{};

//---------------------------------------------------------------
// residual has to be same type as state for now

template<typename residual_type, typename enable = void>
struct isLegitimateImplicitResidualType
  : isLegitimateImplicitStateType<residual_type>{};

  
//---------------------------------------------------------------
// (for now) only enable IFF Jacobian is our matrix class wrapper

template<typename jacobian_type, typename enable = void>
struct isLegitimateJacobianType : std::false_type{};

template<typename jacobian_type>
struct isLegitimateJacobianType<jacobian_type,
       typename std::enable_if<
	 core::meta::is_coreMatrixWrapper<jacobian_type>::value
	 >::type
       > : std::true_type{};


//---------------------------------------------------------------
  

} // namespace meta
} // namespace core
#endif
