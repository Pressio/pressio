
#ifndef ODE_META_HPP_
#define ODE_META_HPP_

#include <type_traits>
#include "meta/core_meta_basic.hpp"
#include "vector/core_vector_meta.hpp"
#include "vector/core_vector_traits.hpp"
#include "matrix/core_matrix_meta.hpp"

namespace ode{
namespace meta {


template<typename time_type, typename enable = void>
struct isLegitimateTimeType : std::false_type{};

template<typename time_type>
struct isLegitimateTimeType<
  time_type,
  typename std::enable_if<std::is_floating_point<time_type>::value>::type
  > : std::true_type{};

//---------------------------------------------------------------
//---------------------------------------------------------------

// For explicit methods, things are easier so the admissible types
// for states is kind of anything as long as has [] operator.
// Because we just wrap it inside, without needing to do much else.
template<typename state_type, typename enable = void>
struct isLegitimateExplicitStateType : std::false_type{};

template<typename state_type>
struct isLegitimateExplicitStateType<state_type,
  typename std::enable_if<core::details::traits<state_type>::isVector==1>::type
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

// template<typename residual_type>
// struct isLegitimateResidualType<
//   residual_type,
//   typename std::enable_if<!std::is_void<residual_type>::value &&
// 			  !std::is_integral<residual_type>::value &&
// 			  !std::is_floating_point<residual_type>::value
// 			  >::type
//   > : std::true_type{};
//---------------------------------------------------------------


// For implicit methods, things are more complicated. We need to solve
// linear or non linear systems. So allow only if the state is wrapped
// with out core::vector class. 
template<typename state_type, typename enable = void>
struct isLegitimateImplicitStateType : std::false_type{};

template<typename state_type>
struct isLegitimateImplicitStateType<state_type,
				     typename std::enable_if<
				       state_type::isVector==1>::type
				     > : std::true_type{};

//---------------------------------------------------------------

// residual has to be same type as state for now
template<typename residual_type, typename enable = void>
struct isLegitimateImplicitResidualType
  : isLegitimateImplicitStateType<residual_type>{};

  
//---------------------------------------------------------------

// (for now) only enable IFF it is our matrix class 
template<typename jacobian_type, typename enable = void>
struct isLegitimateJacobianType : std::false_type{};

template<typename jacobian_type>
struct isLegitimateJacobianType<jacobian_type,
				typename std::enable_if<jacobian_type::isMatrix==1 &&
							 jacobian_type::isEigen==1 &&
							 jacobian_type::isDense==1
							>::type
				> : std::true_type{};

template<typename jacobian_type>
struct isLegitimateJacobianType<jacobian_type,
				typename std::enable_if<jacobian_type::isMatrix==1 &&
							jacobian_type::isEigen==1 &&
							jacobian_type::isSparse==1
							>::type
				> : std::true_type{};


//---------------------------------------------------------------
  

} // namespace meta
} // namespace core
#endif
