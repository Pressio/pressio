
#ifndef ODE_META_EXPLICIT_HPP_
#define ODE_META_EXPLICIT_HPP_

#include "meta/core_meta_basic.hpp"
#include "vector/core_vector_traits.hpp"

namespace ode{
namespace meta {
  
// For explicit methods, things are easy so the admissible types
// for STATE is almost anything as long as its has [] operator.
// Because we just wrap it inside, without needing to do much else.
template<typename state_type, typename enable = void>
struct isLegitimateExplicitStateType : std::false_type{};

template<typename state_type>
struct isLegitimateExplicitStateType<state_type,
 typename std::enable_if<
   core::meta::is_coreVector<state_type>::value
   >::type
  > : std::true_type{};

template<typename state_type>
struct isLegitimateExplicitStateType<state_type,
  typename std::enable_if<
    core::meta::is_vectorStdLib<state_type>::value
    >::type
  > : std::true_type{};

template<typename state_type>
struct isLegitimateExplicitStateType<state_type,
  typename std::enable_if<
    core::meta::is_vectorEigen<state_type>::value
    >::type
  > : std::true_type{};

template<typename state_type>
struct isLegitimateExplicitStateType<state_type,
  typename std::enable_if<
    core::meta::is_vectorEpetra<state_type>::value
    >::type
  > : std::true_type{};
  
//---------------------------------------------------------------
// residual has to be same type as state for now
template<typename residual_type, typename enable = void>
struct isLegitimateExplicitResidualType
  : isLegitimateExplicitStateType<residual_type>{};


} // namespace meta
} // namespace core
#endif
