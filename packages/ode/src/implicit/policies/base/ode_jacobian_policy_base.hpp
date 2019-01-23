
#ifndef ODE_POLICIES_BASE_JACOBIAN_POLICY_BASE_HPP_
#define ODE_POLICIES_BASE_JACOBIAN_POLICY_BASE_HPP_

#include "../../../ode_ConfigDefs.hpp"

namespace rompp{ namespace ode{ namespace policy{

template <typename derived_t>
struct JacobianPolicyBase
  : private core::details::CrtpBase<
	JacobianPolicyBase<derived_t>>{

  using this_t = JacobianPolicyBase<derived_t>;

  /* for now empty, because policies overload (),
   * but we can fill in other methods if needed
   * that all derived policies need to implement
   */

private:
  friend derived_t;
  friend core::details::CrtpBase<this_t>;

  JacobianPolicyBase() = default;
  ~JacobianPolicyBase() = default;

};//end class

}}}//end namespace rompp::ode::policy
#endif
