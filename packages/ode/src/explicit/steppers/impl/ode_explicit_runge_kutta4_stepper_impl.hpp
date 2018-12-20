
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_RUNGEKUTTA4_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_RUNGEKUTTA4_STEPPER_IMPL_HPP_

#include "../base/ode_explicit_stepper_base.hpp"

namespace rompp{ namespace ode{ namespace impl{

template<typename state_type,
	 typename model_type,
	 typename ode_residual_type,
	 typename residual_policy_type
	 >
class ExplicitRungeKutta4StepperImpl<state_type,
				     model_type,
				     ode_residual_type,
				     residual_policy_type>
  : public ExplicitStepperBase<
  ExplicitRungeKutta4StepperImpl<state_type,
				 model_type,
				 ode_residual_type,
				 residual_policy_type> >,
    protected OdeStorage<state_type, ode_residual_type, 1, 4>,
    protected ExpOdeAuxData<model_type, residual_policy_type>
{

  static_assert( meta::is_legitimate_explicit_residual_policy<
		 residual_policy_type>::value ||
		 meta::is_explicit_runge_kutta4_residual_standard_policy<
		 residual_policy_type>::value,
"EXPLICIT RUNGEKUTTA4 RESIDUAL_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");

  using this_t = ExplicitRungeKutta4StepperImpl<
    state_type, model_type, ode_residual_type, residual_policy_type>;
  using stepper_base_t = ExplicitStepperBase<this_t>;
  using storage_base_t = OdeStorage<state_type, ode_residual_type, 1, 4>;
  using auxdata_base_t = ExpOdeAuxData<model_type, residual_policy_type>;
  using scalar_type  = typename core::details::traits<state_type>::scalar_t;
  using scalar_t2  = typename core::details::traits<ode_residual_type>::scalar_t;
  static_assert(std::is_same<scalar_type, scalar_t2>::value,
		"Not maching scalar types");

protected:
  template < typename T1 = model_type,
  	     typename T2 = residual_policy_type,
	     typename T3 = state_type,
	     typename T4 = ode_residual_type,
	     typename... Args>
  ExplicitRungeKutta4StepperImpl(const T1 & model,
				 const T2 & res_policy_obj,
				 const T3 & y0,
				 const T4 & r0,
				 Args&&... rest)
    : storage_base_t(y0, r0), auxdata_base_t(model, res_policy_obj)
  {
    //make sure there is something in what is passed,
    //otherwise the helper data structures are emtpy
    assert( !y0.empty() );
    assert( !r0.empty() );
  }

  template < typename T1 = residual_policy_type,
	     typename T2 = state_type,
	     typename T3 = ode_residual_type>
  ExplicitRungeKutta4StepperImpl(const T1 & res_policy_obj,
				 const T2 & y0,
				 const T3 & r0)
    : storage_base_t(y0, r0), auxdata_base_t(res_policy_obj)
  {
    //make sure there is something in what is passed,
    //otherwise the helper data structures are emtpy
    assert( !y0.empty() );
    assert( !r0.empty() );
  }

  ExplicitRungeKutta4StepperImpl() = delete;
  ~ExplicitRungeKutta4StepperImpl() = default;

public:

  template<typename step_t>
  void operator()(state_type & y, scalar_type t,
		  scalar_type dt, step_t step){

    auto & ytmp = this->auxStates_[0];

    const scalar_type dt_half = dt / static_cast< scalar_type >(2);
    const scalar_type t_phalf = t + dt_half;
    const scalar_type dt6 = dt / static_cast< scalar_type >( 6 );
    const scalar_type dt3 = dt / static_cast< scalar_type >( 3 );

    // ----------
    // stage 1:
    // ----------
    (*this->residual_obj_)(y, this->auxRHS_[0], *this->model_, t);
    ytmp = y + this->auxRHS_[0]*dt_half;
    //ytmp.template inPlaceOp<add_op_t>(1.0,  y, dt_half, this->auxRHS_[0]);

    // ----------
    // stage 2:
    // ----------
    (*this->residual_obj_)(ytmp, this->auxRHS_[1], *this->model_, t_phalf);
    ytmp = y + this->auxRHS_[1]*dt_half;
    //ytmp.template inPlaceOp<add_op_t>(1.0, y, dt_half, this->auxRHS_[1]);

    // ----------
    // stage 3:
    // ----------
    (*this->residual_obj_)(ytmp, this->auxRHS_[2], *this->model_, t_phalf);
    ytmp = y + this->auxRHS_[2]*dt;
    //ytmp.template inPlaceOp<add_op_t>(1.0, y, dt, this->auxRHS_[2]);

    // ----------
    // stage 4:
    // ----------
    (*this->residual_obj_)(ytmp, this->auxRHS_[3], *this->model_, t + dt);
    //y_n += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    y += dt6*this->auxRHS_[0] + dt3*this->auxRHS_[1] + dt3*this->auxRHS_[2] + dt6*this->auxRHS_[3];
    // y.template inPlaceOp<add_op_t>(1.0, dt6, auxRHS_[0],
    // 				   dt3, auxRHS_[1],
    // 				   dt3, auxRHS_[2],
    // 				   dt6, auxRHS_[3]);
  }//end doStep

private:
  friend stepper_base_t;

}; //end class

}}}//end namespace rompp::ode::impl
#endif
