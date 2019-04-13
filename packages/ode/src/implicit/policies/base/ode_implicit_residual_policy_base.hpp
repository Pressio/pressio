
#ifndef ODE_POLICIES_BASE_IMPLICIT_RESIDUAL_POLICY_BASE_HPP_
#define ODE_POLICIES_BASE_IMPLICIT_RESIDUAL_POLICY_BASE_HPP_

#include "../../../ode_ConfigDefs.hpp"

namespace rompp{ namespace ode{ namespace policy{

template <typename derived_t>
struct ImplicitResidualPolicyBase
  : private core::details::CrtpBase<
  ImplicitResidualPolicyBase<derived_t>>{

  using this_t = ImplicitResidualPolicyBase<derived_t>;

  /* for now empty, because policies overload (),
   * but we can fill in other methods if needed
   * that all derived policies need to implement
   */

private:
  /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  template<typename DummyType> struct dummy{using type = DummyType;};
  friend typename dummy<derived_t>::type;

  friend core::details::CrtpBase<this_t>;

  ImplicitResidualPolicyBase() = default;
  ~ImplicitResidualPolicyBase() = default;

};//end class


}}}//end namespace rompp::ode::policy
#endif



//-----------------------------------------------------
// when computing TIME residudal but we also
// need  previous RHS not just previous states
//-----------------------------------------------------
// template <typename state_type,
// 	    typename residual_type,
// 	    typename model_type,
// 	    typename scalar_type,
// 	    int T = numAuxRHS,
// 	    typename std::enable_if<T!=0>::type * = nullptr>
// void operator()(const state_type & y,
// 	       residual_type & R,
// 	       const std::array<state_type, numAuxStates> & auxYs,
// 	       const std::array<residual_type, T> & auxRHSs,
// 	       model_type & model,
// 	       scalar_type t,
// 	       scalar_type dt)
// {
//   this->underlying()(y, R, auxYs, auxRHSs, model, t, dt);
// }
//-----------------------------------------------------
