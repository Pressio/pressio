
#ifndef ODE_EXPLICIT_STEPPER_BASE_HPP_
#define ODE_EXPLICIT_STEPPER_BASE_HPP_

#include "ode_ConfigDefs.hpp"
#include "ode_stepper_traits.hpp"


namespace ode{

template<typename stepper_type>
class explicit_stepper_base
{
public:

  using state_t = typename ode::details::traits<stepper_type>::state_t;
  using sc_t = typename ode::details::traits<stepper_type>::scalar_t;
  using order_t = typename ode::details::traits<stepper_type>::order_t; 
  static constexpr order_t order_value = ode::details::traits<stepper_type>::order_value;
  using time_t = ode::details::time_type;
  using resizer_t = typename ode::details::traits<stepper_type>::resizer_t;
  using rhs_t = typename ode::details::traits<stepper_type>::rhs_t;

  // (de)constructors
  explicit_stepper_base(){}
  virtual ~explicit_stepper_base(){}

  // methods
  order_t order() const{  return order_value; }

  template< typename functor_type>
  void doStep( functor_type & functor, state_t &inout, time_t t, time_t dt ){
    this->stepper()->doStepImpl( functor, inout, t, dt );
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


  // template< typename functor_type>
  // void do_step( functor_type & functor, state_type &inout, time_type t, time_type dt ){
  //   functor( inout, rhs_, t );
  //   this->stepper()->do_step_impl( functor, inout, rhs_, inout, t, dt );
  // }


  // typedef State state_type;
  // typedef Value value_type;
  // typedef Deriv deriv_type;
  // typedef Time time_type;
  // typedef Resizer resizer_type;
  // typedef Stepper stepper_type;
  // typedef stepper_tag stepper_category;
  // typedef algebra_stepper_base< Algebra , Operations > algebra_stepper_base_type;
  // typedef typename algebra_stepper_base_type::algebra_type algebra_type;
  // typedef typename algebra_stepper_base_type::operations_type operations_type;
  // typedef unsigned short order_type;
  // typedef state_wrapper< state_type > wrapped_state_type;
  // typedef state_wrapper< deriv_type > wrapped_deriv_type;
