
#ifndef ODE_EXPLICIT_STEPPER_BASE_HPP_
#define ODE_EXPLICIT_STEPPER_BASE_HPP_

#include "ode_ConfigDefs.hpp"

// template<
//   typename stepper_type,
//   unsigned int order,
//   typename state_type,
//   typename deriv_type = state_type,
//   typename scalar_type = typename core::defaultTypes::scalar_t,
//   typename time_type = typename core::defaultTypes::scalar_t
// >

template<typename stepper_type>
class explicit_stepper_base
{
public:

  using state_t = typename ode::details::traits<stepper_type>::state_t;
  using sc_t = typename ode::details::traits<stepper_type>::scalar_t;
  using order_type = typename ode::details::traits<stepper_type>::order_t;
  static constexpr order_type order_value = ode::details::traits<stepper_type>::order_value;

  // (de)constructors
  explicit_stepper_base(){}
  virtual ~explicit_stepper_base(){}

  // methods
  order_type order() const{  return order_value; }

  template< typename functor_type>
  void do_step( functor_type & functor, state_type &inout, time_type t, time_type dt )
  {
    this->stepper()->do_step_impl( functor, inout, inout, t, dt );
  }

private:
    stepper_type * stepper( ){
        return static_cast< stepper_type* >( this );
    }
    const stepper_type * stepper( void ) const{
        return static_cast< const stepper_type* >( this );
    }

};

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
