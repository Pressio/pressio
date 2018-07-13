
#ifndef ODE_EXPLICIT_POLICIES_META_HPP_
#define ODE_EXPLICIT_POLICIES_META_HPP_

#include "../base/ode_explicit_residual_policy_base.hpp"

namespace ode{
namespace meta {

//-----------------------------------------------------------------
// METAF FOR ADMISSIBLE EXPLICIT EULER RESIDUAL
//-----------------------------------------------------------------
template<typename policy_t, typename enable = void>
struct isLegitimateExplicitResidualPolicy
  : std::false_type{};
  
template <typename policy_t>
struct isLegitimateExplicitResidualPolicy<
  policy_t,
  typename std::enable_if<
    core::meta::publiclyInheritsFrom<
      policy_t,
      ode::policy::explicitResidualPolicyBase<policy_t>
      >::value
    >::type
  > : std::true_type{};


} // namespace meta
} // namespace core
#endif
