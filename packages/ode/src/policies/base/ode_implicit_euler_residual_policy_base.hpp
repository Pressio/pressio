
#ifndef ODE_IMPLICIT_EULER_RESIDUAL_POLICY_BASE_HPP_
#define ODE_IMPLICIT_EULER_RESIDUAL_POLICY_BASE_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace policy{
  
template<typename derived_type,
	 typename state_type,
	 typename residual_type,
	 typename model_type,
	 typename time_type>
class implicitEulerResidualPolicyBase{
public:
  void compute(const state_type & y, const state_type & ynm1,
	       residual_type & R, model_type & model,
	       time_type t, time_type dt){
    this->underlying()->computeImpl(y,ynm1,R,model,t,dt);
  } 

private:
  friend derived_type; 
  implicitEulerResidualPolicyBase() = default;
  ~implicitEulerResidualPolicyBase() = default;
  
  derived_type & underlying(){
    return static_cast<derived_type& >( *this );
  }
  derived_type const & underlying() const{
    return static_cast<derived_type const & >( *this );
  }
};

}//end namespace polices
}//end namespace ode  
#endif 
