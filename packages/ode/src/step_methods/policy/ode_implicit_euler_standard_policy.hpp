
#ifndef ODE_IMPLICIT_EULER_STANDARD_POLICY_HPP_
#define ODE_IMPLICIT_EULER_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "../impl/ode_euler_implicit_impl.hpp"
#include "ode_implicit_euler_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type,
	 typename residual_type,
	 typename model_type,
	 typename time_type>
class implicitEulerStandardResidual
  : public implicitEulerResidualPolicyBase<implicitEulerStandardResidual<state_type,
									  residual_type,
									  model_type,
									  time_type>,
					   state_type, residual_type,
					   model_type, time_type>
{
public:
  void compute(const state_type & y, const state_type & ynm1,
	       residual_type & R, model_type & model,
	       time_type t, time_type dt)
  {
    // first eval RHS
    model.residual(y,R,t);
    // then fix residual based on time stepping features
    ode::impl::implicit_euler_residual_impl(y, ynm1, R, dt);
  }
};


template<typename state_type,
	 typename jacobian_type,
	 typename model_type, 
	 typename time_type>
class implicitEulerStandardJacobian
  : public implicitEulerJacobianPolicyBase<implicitEulerStandardJacobian<state_type,
									 jacobian_type,
									 model_type,
									 time_type>,
					   state_type, jacobian_type,
					   model_type,time_type
					   > 
{
public:
  void compute(const state_type & y, jacobian_type & J,
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

