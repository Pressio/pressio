
#ifndef ODE_META_META_IMPLICIT_HPP_
#define ODE_META_META_IMPLICIT_HPP_

#include "meta/core_meta_basic.hpp"
#include "vector/core_vector_meta.hpp"
#include "matrix/core_matrix_meta.hpp"

namespace ode{
namespace meta {

// For implicit methods, things are more complicated. We need to solve
// linear or non linear systems. So allow ONLY if the state is wrapped
// with out core::Vector<> class. 

template<typename state_type, typename enable = void>
struct isLegitimateImplicitStateType : std::false_type{};

template<typename state_type>
struct isLegitimateImplicitStateType<state_type,
        typename std::enable_if<
	  core::meta::is_core_vector_wrapper<state_type>::value
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
	 core::meta::is_core_matrix_wrapper<jacobian_type>::value
	 >::type
       > : std::true_type{};
  

} // namespace meta
} // namespace core
#endif
