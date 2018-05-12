
#ifndef TIMEINTEGRATOR_EXPLICITSTEPPERBASE_HPP
#define TIMEINTEGRATOR_EXPLICITSTEPPERBASE_HPP

#include "timeIntegrator_ConfigDefs.hpp"

template<
  class stepper_type,
  unsigned int order_type,
  class state_type,
  class deriv_type = state_type,
  class scalar_type = typename timeIntegrator::details::defaultTypes::scalar_t,
  class time_type = typename timeIntegrator::details::defaultTypes::scalar_t
>
class explicit_stepper_base
{
public:
  static constexpr order_type order_value = Order;

  // (de)constructors
  explicit_stepper_base(){}
  virtual ~explicit_stepper_base(){}

  // methods
  order_type order() const{  return order_value; }

  template< class system_type>
  void do_step( system_type & functor,
		state_type &inout,
		time_type t,
		time_type dt )
  {
    functor( inout, rhs_, t );
    this->stepper()->step_impl( functor, inout, rhs_, inout, t, dt );
  }

private:
    stepper_type * stepper( ){
        return static_cast< stepper_type* >( this );
    }
    const stepper_type * stepper( void ) const{
        return static_cast< const stepper_type* >( this );
    }

protected:
    deriv_type rhs_;

};



  // typedef explicit_stepper_base< Stepper, Order , State , Value , Deriv , Time , Algebra , Operations , Resizer > internal_stepper_base_type;

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
