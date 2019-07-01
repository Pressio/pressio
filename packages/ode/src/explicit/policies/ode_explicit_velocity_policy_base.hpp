
#ifndef ODE_POLICIES_BASE_EXPLICIT_VELOCITY_POLICY_BASE_HPP_
#define ODE_POLICIES_BASE_EXPLICIT_VELOCITY_POLICY_BASE_HPP_

#include "../../ode_fwd.hpp"

namespace rompp{ namespace ode{ namespace policy{

template <typename derived_t>
class ExplicitVelocityPolicyBase
{

  using this_t = ExplicitVelocityPolicyBase<derived_t>;

  /* for now empty, because policies overload (),
   * but we can fill in other methods if needed
   * that all derived policies need to implement
   */

public:
  ExplicitVelocityPolicyBase() = default;
  ~ExplicitVelocityPolicyBase() = default;

};//end class

}}}//end namespace rompp::ode::policy
#endif
