
#ifndef ODE_POLICIES_BASE_JACOBIAN_POLICY_BASE_HPP_
#define ODE_POLICIES_BASE_JACOBIAN_POLICY_BASE_HPP_

#include "../../../ode_ConfigDefs.hpp"

namespace pressio{ namespace ode{ namespace policy{

template <typename derived_t>
struct JacobianPolicyBase
  : private utils::details::CrtpBase<
	JacobianPolicyBase<derived_t>>{

  using this_t = JacobianPolicyBase<derived_t>;

  /* for now empty, because policies overload (),
   * but we can fill in other methods if needed
   * that all derived policies need to implement
   */

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_t>::type;

  friend utils::details::CrtpBase<this_t>;

  JacobianPolicyBase() = default;
  ~JacobianPolicyBase() = default;

};//end class

}}}//end namespace pressio::ode::policy
#endif
