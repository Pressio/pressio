
#ifndef ODE_EXPLICIT_POLICIES_IS_LEGITIMATE_EXPLICIT_RESID_POL_HPP_
#define ODE_EXPLICIT_POLICIES_IS_LEGITIMATE_EXPLICIT_RESID_POL_HPP_

#include "ode_explicit_residual_standard_policy.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<typename policy_t, typename enable = void>
struct is_legitimate_explicit_residual_policy
  : std::false_type{};

template <typename policy_t>
struct is_legitimate_explicit_residual_policy<
  policy_t,
  typename std::enable_if<
    ::rompp::mpl::publicly_inherits_from<
      policy_t,
      ode::policy::ExplicitResidualPolicyBase<policy_t>
      >::value
    >::type
  > : std::true_type{};


}}}//end namespace rompp::core::meta
#endif
