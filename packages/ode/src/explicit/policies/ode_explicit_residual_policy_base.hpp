
#ifndef ODE_POLICIES_BASE_EXPLICIT_RESIDUAL_POLICY_BASE_HPP_
#define ODE_POLICIES_BASE_EXPLICIT_RESIDUAL_POLICY_BASE_HPP_

#include "../../ode_forward_declarations.hpp"

namespace rompp{ namespace ode{ namespace policy{

template <typename derived_t>
class ExplicitResidualPolicyBase
  : private core::details::CrtpBase<
  ExplicitResidualPolicyBase<derived_t>>{

  using this_t = ExplicitResidualPolicyBase<derived_t>;

  /* for now empty, because policies overload (),
   * but we can fill in other methods if needed
   * that all derived policies need to implement
   */

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_t>::type;

  friend core::details::CrtpBase<this_t>;

  ExplicitResidualPolicyBase() = default;
  ~ExplicitResidualPolicyBase() = default;

};//end class

}}}//end namespace rompp::ode::policy
#endif
