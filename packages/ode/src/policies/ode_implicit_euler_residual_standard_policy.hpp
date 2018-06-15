
#ifndef ODE_IMPLICIT_EULER_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_IMPLICIT_EULER_RESIDUAL_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../step_methods/impl/ode_euler_implicit_residual_impl.hpp"
#include "./base/ode_implicit_euler_residual_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type, typename residual_type,
	 typename model_type, typename time_type>
class implicitEulerStandardResidual
  : public implicitEulerResidualPolicyBase<
             implicitEulerStandardResidual,state_type,
                 residual_type,model_type, time_type>
{
  implicitEulerStandardResidual() = default;
  ~implicitEulerStandardResidual() = default;  
private:
  void computeImpl(const state_type & y, const state_type & ynm1,
		   residual_type & R, model_type & model,
		   time_type t, time_type dt)
  {
    // first eval RHS
    model.residual(y,R,t);
    // then fix residual based on time stepping features
    ode::impl::implicit_euler_residual_impl(y, ynm1, R, dt);
  }
private:
  friend implicitEulerResidualPolicyBase<implicitEulerStandardResidual,
					 state_type,residual_type,
					 model_type, time_type>;
};

  
}//end namespace polices
}//end namespace ode  
#endif 

