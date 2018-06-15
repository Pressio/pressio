
#ifndef ODE_IMPLICIT_EULER_JACOBIAN_STANDARD_POLICY_HPP_
#define ODE_IMPLICIT_EULER_JACOBIAN_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "./step_methods/impl/ode_euler_implicit_jacobian_impl.hpp"
#include "./base/ode_implicit_euler_jacobian_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type,
	 typename jacobian_type,
	 typename model_type, 
	 typename time_type>
class implicitEulerStandardJacobian
  : public implicitEulerJacobianPolicyBase<
  implicitEulerStandardJacobian,state_type,
  jacobian_type,model_type,time_type>
{
  implicitEulerStandardJacobian() = default;
  ~implicitEulerStandardJacobian() = default;
private:
  void computeImpl(const state_type & y, jacobian_type & J,
		   model_type & model, time_type t, time_type dt){
    // first eval jac
    model.jacobian(y,J,t);
    // then fix it based on time stepping features
    ode::impl::implicit_euler_jacobian_impl(y, J, dt);
  }
};

  
}//end namespace polices
}//end namespace ode  
#endif 
