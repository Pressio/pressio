
#ifndef ODE_IMPLICIT_EULER_POLICY_BASE_HPP_
#define ODE_IMPLICIT_EULER_POLICY_BASE_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace policy{
  
template<typename derived_type,
	 typename state_type, typename residual_type,
	 typename model_type, typename time_type>
class implicitEulerResidualPolicyBase{
public:
  derived_type * underlying(){
    return static_cast< derived_type* >( this );
  }
  const derived_type * underlying() const{
    return static_cast< const derived_type* >( this );
  }

  void compute(const state_type & y, const state_type & ynm1,
	       residual_type & R, model_type & model,
	       time_type t, time_type dt){
    this->underlying()->computeImpl(y,ynm1,R,model,t,dt);
  } 
};



template <typename derived_type,
	  typename state_type,
	  typename jacobian_type,
	  typename model_type,
	  typename time_type>
class implicitEulerJacobianPolicyBase{
public:
  derived_type * underlying(){
    return static_cast< derived_type* >( this );
  }
  const derived_type * underlying() const{
    return static_cast< const derived_type* >( this );
  }

  void compute(const state_type & y, jacobian_type & J,
	       model_type & model, time_type t, time_type dt){
    this->underlying()->computeImpl(y,J,model,t,dt);
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
