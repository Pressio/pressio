
#ifndef ODE_EXPLICIT_RUNGEKUTTA4_POLICY_BASE_HPP_
#define ODE_EXPLICIT_RUNGEKUTTA4_POLICY_BASE_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace policy{
    
template <typename derived_type,
	  typename state_type, typename residual_type,
	  typename model_type, typename time_type
	  >
class explicitRungeKutta4ResidualPolicyBase{
public:
  derived_type * underlying(){
    return static_cast< derived_type* >( this );
  }
  const derived_type * underlying() const{
    return static_cast< const derived_type* >( this );
  }
  
  void compute(const state_type & y, residual_type & R,
	       model_type & model, time_type t){
    this->underlying()->computeImpl(y, R, model, t);
  } 
};
  
}//end namespace polices
}//end namespace ode  
#endif 
