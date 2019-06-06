
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_RUNGEKUTTA4_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_RUNGEKUTTA4_STEPPER_IMPL_HPP_

#include "../ode_explicit_stepper_base.hpp"

namespace rompp{ namespace ode{ namespace impl{

template<
  typename scalar_type,
  typename state_type,
  typename model_type,
  typename ode_residual_type,
  typename residual_policy_type,
  typename ops_t
  >
class ExplicitRungeKutta4StepperImpl<scalar_type,
				     state_type,
				     model_type,
				     ode_residual_type,
				     residual_policy_type,
				     ops_t>
  : protected OdeStorage<state_type, ode_residual_type, 1, 4>,
    protected ExplicitOdeAuxData<model_type, residual_policy_type>
{

  static_assert( meta::is_legitimate_explicit_residual_policy<
		 residual_policy_type>::value ||
		 meta::is_explicit_runge_kutta4_residual_standard_policy<
		 residual_policy_type>::value,
"EXPLICIT RUNGEKUTTA4 RESIDUAL_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");

  using this_t = ExplicitRungeKutta4StepperImpl< scalar_type,
						 state_type, model_type,
						 ode_residual_type,
						 residual_policy_type,
						 ops_t>;

  using storage_base_t = OdeStorage<state_type, ode_residual_type, 1, 4>;
  using auxdata_base_t = ExplicitOdeAuxData<model_type, residual_policy_type>;

  using storage_base_t::auxRHS_;
  using storage_base_t::auxStates_;
  using auxdata_base_t::model_;
  using auxdata_base_t::residual_obj_;

public:
  ExplicitRungeKutta4StepperImpl(const model_type & model,
				 const residual_policy_type & res_policy_obj,
				 const state_type & y0,
				 const ode_residual_type & r0)
    : storage_base_t(y0, r0), auxdata_base_t(model, res_policy_obj){}

  ExplicitRungeKutta4StepperImpl(const residual_policy_type & res_policy_obj,
				 const state_type & y0,
				 const ode_residual_type & r0)
    : storage_base_t(y0, r0), auxdata_base_t(res_policy_obj){}

  ExplicitRungeKutta4StepperImpl() = delete;
  ~ExplicitRungeKutta4StepperImpl() = default;

public:

  // if user does NOT provide ops, then user wrapper ops
  template<typename step_t,
  	   typename T = ops_t,
  	   mpl::enable_if_t<
  	     std::is_void<T>::value
  	     > * = nullptr
  	   >
  void doStep(state_type & y,
  	      scalar_type t,
  	      scalar_type dt,
  	      step_t step){

    auto & ytmp = this->auxStates_[0];

    const scalar_type dt_half = dt / static_cast< scalar_type >(2);
    const scalar_type t_phalf = t + dt_half;
    const scalar_type dt6 = dt / static_cast< scalar_type >( 6 );
    const scalar_type dt3 = dt / static_cast< scalar_type >( 3 );

    constexpr auto one  = ::rompp::core::constants::one<scalar_type>();

    // stage 1: ytmp = y + auxRHS_[0]*dt_half;
    (*this->residual_obj_)(y, this->auxRHS_[0], *this->model_, t);
    ::rompp::core::ops::do_update(ytmp, y, one, this->auxRHS_[0], dt_half);

    // stage 2: ytmp = y + auxRHS_[1]*dt_half;
    (*this->residual_obj_)(ytmp, this->auxRHS_[1], *this->model_, t_phalf);
    ::rompp::core::ops::do_update(ytmp, y, one, this->auxRHS_[1], dt_half);

    // stage 3: ytmp = y + auxRHS_[2]*dt;
    (*this->residual_obj_)(ytmp, this->auxRHS_[2], *this->model_, t_phalf);
    ::rompp::core::ops::do_update(ytmp, y, one, this->auxRHS_[2], dt);

    // stage 4: y_n += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    (*this->residual_obj_)(ytmp, this->auxRHS_[3], *this->model_, t + dt);
    ::rompp::core::ops::do_update(y, one,
  				  this->auxRHS_[0], dt6,
  				  this->auxRHS_[1], dt3,
  				  this->auxRHS_[2], dt3,
  				  this->auxRHS_[3], dt6);
  }//end doStep


  // if user provides ops, then use them
  template<typename step_t,
	   typename T = ops_t,
	   mpl::enable_if_t<
	     std::is_void<T>::value == false
	     > * = nullptr
	   >
  void doStep(state_type & y,
	      scalar_type t,
	      scalar_type dt,
	      step_t step){
    using op = typename ops_t::update_op;

    auto & ytmp = auxStates_[0];

    const scalar_type dt_half = dt / static_cast< scalar_type >(2);
    const scalar_type t_phalf = t + dt_half;
    const scalar_type dt6 = dt / static_cast< scalar_type >( 6 );
    const scalar_type dt3 = dt / static_cast< scalar_type >( 3 );

    constexpr auto one  = ::rompp::core::constants::one<scalar_type>();

    // stage 1: ytmp = y + auxRHS_[0]*dt_half;
    (*residual_obj_)(y, auxRHS_[0], *model_, t);
    op::do_update(*ytmp.data(), *y.data(), one, *auxRHS_[0].data(), dt_half);

    // stage 2: ytmp = y + auxRHS_[1]*dt_half;
    (*residual_obj_)(ytmp, auxRHS_[1], *model_, t_phalf);
    op::do_update(*ytmp.data(), *y.data(), one, *auxRHS_[1].data(), dt_half);

    // stage 3: ytmp = y + auxRHS_[2]*dt;
    (*residual_obj_)(ytmp, auxRHS_[2], *model_, t_phalf);
    op::do_update(*ytmp.data(), *y.data(), one, *auxRHS_[2].data(), dt);

    // stage 4: y_n += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    (*residual_obj_)(ytmp, auxRHS_[3], *model_, t + dt);
    op::do_update(*y.data(), one,
		  *auxRHS_[0].data(), dt6,
		  *auxRHS_[1].data(), dt3,
		  *auxRHS_[2].data(), dt3,
		  *auxRHS_[3].data(), dt6);
  }//end doStep

}; //end class

}}}//end namespace rompp::ode::impl
#endif
