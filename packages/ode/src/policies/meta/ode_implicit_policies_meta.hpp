
#ifndef ODE_IMPLICIT_POLICIES_META_HPP_
#define ODE_IMPLICIT_POLICIES_META_HPP_

#include "../base/ode_implicit_residual_policy_base.hpp"
#include "../base/ode_jacobian_policy_base.hpp"

namespace ode{
namespace meta {

//-----------------------------------------------------------------
// METAF FOR ADMISSIBLE JACOBIAN POLICY
//-----------------------------------------------------------------
template<typename policy_t, typename enable = void>
struct isLegitimateImplicitJacobianPolicy : std::false_type{};
  
template <typename policy_t>
struct isLegitimateImplicitJacobianPolicy<
  policy_t,
  typename std::enable_if<
    core::meta::publiclyInheritsFrom<
      policy_t,
      ode::policy::jacobianPolicyBase<policy_t>
      >::value
    >::type
  > : std::true_type{};
//-----------------------------------------------------------------
    
} // namespace meta
} // namespace core
#endif
