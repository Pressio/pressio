
#ifndef ODE_POLICIES_META_EXPLICIT_POLICIES_META_HPP_
#define ODE_POLICIES_META_EXPLICIT_POLICIES_META_HPP_

#include "../base/ode_explicit_residual_policy_base.hpp"

namespace ode{
namespace meta {

//-----------------------------------------------------------------
// METAF FOR ADMISSIBLE EXPLICIT EULER RESIDUAL
//-----------------------------------------------------------------
template<typename policy_t, typename enable = void>
struct is_legitimate_explicit_residual_policy
  : std::false_type{};
  
template <typename policy_t>
struct is_legitimate_explicit_residual_policy<
  policy_t,
  typename std::enable_if<
    core::meta::publicly_inherits_from<
      policy_t,
      ode::policy::ExplicitResidualPolicyBase<policy_t>
      >::value
    >::type
  > : std::true_type{};


} // namespace meta
} // namespace core
#endif
