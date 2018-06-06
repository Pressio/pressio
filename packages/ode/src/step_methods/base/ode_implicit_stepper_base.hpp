
#ifndef ODE_IMPLICIT_STEPPER_BASE_HPP_
#define ODE_IMPLICIT_STEPPER_BASE_HPP_

#include "ode_ConfigDefs.hpp"
#include "../ode_stepper_traits.hpp"


namespace ode{

template<typename stepper_type>
class implicitStepperBase
{
public:
  using state_t = typename ode::details::traits<stepper_type>::state_t;
  using sc_t = typename ode::details::traits<stepper_type>::scalar_t;
  using order_t = typename ode::details::traits<stepper_type>::order_t; 
  using time_t = ode::details::time_type;
  using resizer_t = typename ode::details::traits<stepper_type>::resizer_t;
  using rhs_t = typename ode::details::traits<stepper_type>::rhs_t;
  using jacobian_t = typename ode::details::traits<stepper_type>::jacobian_t;

  static constexpr order_t order_value = ode::details::traits<stepper_type>::order_value;
  
  // (de)constructors
  implicitStepperBase(){}
  ~implicitStepperBase(){}

  // methods
  order_t order() const{  return order_value; }

  void doStep(state_t &inout, time_t t, time_t dt){
    this->stepper()->doStepImpl( inout, t, dt );
  }

  void residual(const state_t & y, state_t & R){
    this->stepper()->residualImpl(y, R);
  }

  void jacobian(const state_t & y, jacobian_t & J){
    this->stepper()->jacobianImpl(y, J);
  }
  
  
private:  
  stepper_type * stepper( ){
    return static_cast< stepper_type* >( this );
  }
  const stepper_type * stepper( void ) const{
    return static_cast< const stepper_type* >( this );
  }

protected:
  resizer_t myResizer_;

};


}//end namespace  
#endif
