
#ifndef ODE_POLICIES_STANDARD_IMPLICIT_JACOBIAN_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_JACOBIAN_STANDARD_POLICY_HPP_

#include "../../../ode_forward_declarations.hpp"
#include "../base/ode_jacobian_policy_base.hpp"
#include "../../ode_jacobian_impl.hpp"
#include "../../../ode_basic_meta.hpp"

namespace rompp{ namespace ode{ namespace policy{

template<typename state_type,
	 typename model_type,
	 typename jacobian_type>
class ImplicitJacobianStandardPolicy<
  state_type, model_type, jacobian_type,
  core::meta::enable_if_t<
    ::rompp::ode::meta::is_legitimate_implicit_state_type<state_type>::value and
    ::rompp::ode::meta::is_legitimate_jacobian_type<jacobian_type>::value
    >
  > : public JacobianPolicyBase<ImplicitJacobianStandardPolicy<
    state_type, model_type, jacobian_type> >{

  using this_t = ImplicitJacobianStandardPolicy<state_type, model_type, jacobian_type>;

public:
  ImplicitJacobianStandardPolicy() = default;
  ~ImplicitJacobianStandardPolicy() = default;

public:
  template <ode::ImplicitEnum method, typename scalar_t>
  void operator()(const state_type & y,
		  jacobian_type & J,
		  const model_type & model,
		  scalar_t t,
		  scalar_t dt)const
  {
    model.jacobian( *y.data(), *J.data(), t);
    ode::impl::time_discrete_jacobian<method>(J, dt);
  }

  template <ode::ImplicitEnum method, typename scalar_t>
  jacobian_type operator()(const state_type & y,
  			   const model_type & model,
  			   scalar_t t,
  			   scalar_t dt)const
  {
    auto nJJ = model.jacobian(*y.data(), t);
    jacobian_type JJ(nJJ);
    ode::impl::time_discrete_jacobian<method>(JJ, dt);
    return JJ;
  }

private:
  friend JacobianPolicyBase<this_t>;

};//end class

}}}//end namespace rompp::ode::policy
#endif
