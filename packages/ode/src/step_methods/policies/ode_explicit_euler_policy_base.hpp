
#ifndef ODE_EXPLICIT_EULER_POLICY_BASE_HPP_
#define ODE_EXPLICIT_EULER_POLICY_BASE_HPP_

#include "ode_ConfigDefs.hpp"


namespace ode{  
namespace policies{  

    
template <typename derived_type>
class explicitEulerResidualPolicyBase : core::crtpBase<derived_type,
						       explicitEulerResidualPolicyBase>
{
public:
  explicitEulerResidualPolicyBase(){}

  template<typename state_type,
	   typename residual_type,
	   typename model_type,
	   typename time_type>
  void compute(const state_type & y, residual_type & R,
	       model_type & model, time_type t){
    this->underlying()->computeImpl(y, R, model, t);
  } 
};

  
}//end namespace polices
}//end namespace ode  
#endif 




// private:
//   derived_type * policy(){
//     return static_cast< derived_type* >( this );
//   }
//   const derived_type * policy() const{
//     return static_cast< const derived_type* >( this );
//   }

// template <typename model_type>
// class policyBase{
// public:
//   policyBase(model_type * model) : model_(model){}  
// protected:
//    model_type * model_;
// };
