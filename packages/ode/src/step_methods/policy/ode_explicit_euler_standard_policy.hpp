#ifndef ODE_EXPLICIT_EULER_STANDARD_POLICY_HPP_
#define ODE_EXPLICIT_EULER_STANDARD_POLICY_HPP_

#include "ode_ConfigDefs.hpp"
#include "ode_explicit_euler_policy_base.hpp"

namespace ode{
namespace policy{

template<typename state_type, typename residual_type,
	 typename model_type, typename time_type>
class explicitEulerStandardResidual
  : public explicitEulerResidualPolicyBase< explicitEulerStandardResidual<state_type,residual_type,
									  model_type,time_type>,
					    state_type, residual_type,model_type, time_type>
{
 public:
  void computeImpl(const state_type & y, residual_type & R,
		   model_type & model, time_type t){
    model.residual(y, R, t);
  }
};
  
}//end namespace polices
}//end namespace ode  
#endif 
