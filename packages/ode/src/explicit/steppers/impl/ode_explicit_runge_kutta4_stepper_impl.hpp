
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_RUNGEKUTTA4_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_RUNGEKUTTA4_STEPPER_IMPL_HPP_

#include "../ode_explicit_stepper_base.hpp"

namespace pressio{ namespace ode{ namespace impl{

template<
  typename scalar_type,
  typename state_type,
  typename system_type,
  typename ode_residual_type,
  typename residual_policy_type,
  typename ops_t
  >
class ExplicitRungeKutta4StepperImpl<scalar_type,
				     state_type,
				     system_type,
				     ode_residual_type,
				     residual_policy_type,
				     ops_t>
{

  static_assert( meta::is_legitimate_explicit_velocity_policy<
		 residual_policy_type>::value ||
		 meta::is_explicit_runge_kutta4_velocity_standard_policy<
		 residual_policy_type>::value,
"EXPLICIT RUNGEKUTTA4 VELOCITY_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");

  using this_t = ExplicitRungeKutta4StepperImpl< scalar_type,
						 state_type, system_type,
						 ode_residual_type,
						 residual_policy_type,
						 ops_t>;

  using state_storage_t = OdeStorage<state_type, 1>;
  using resid_storage_t = OdeStorage<ode_residual_type, 4>;
  using system_wrapper_t = OdeSystemWrapper<system_type>;

  state_storage_t stateAuxStorage_;
  resid_storage_t residAuxStorage_;
  system_wrapper_t sys_;
  const residual_policy_type & policy_;

public:
  ExplicitRungeKutta4StepperImpl(const system_type & model,
  				 const residual_policy_type & res_policy_obj,
  				 const state_type & y0,
  				 const ode_residual_type & r0)
    : stateAuxStorage_{y0}, residAuxStorage_{r0},
      sys_{model}, policy_{res_policy_obj}{}

  ExplicitRungeKutta4StepperImpl() = delete;
  ~ExplicitRungeKutta4StepperImpl() = default;

public:

  // if user does NOT provide ops, then user wrapper ops
  template<
    typename step_t,
    typename T = ops_t,
    typename _state_type = state_type,
    mpl::enable_if_t<
      std::is_void<T>::value
      > * = nullptr
  >
  void doStep(_state_type & y,
  	      scalar_type t,
  	      scalar_type dt,
  	      step_t step){

    auto & ytmp	   = stateAuxStorage_.data_[0];
    auto & auxRhs0 = residAuxStorage_.data_[0];
    auto & auxRhs1 = residAuxStorage_.data_[1];
    auto & auxRhs2 = residAuxStorage_.data_[2];
    auto & auxRhs3 = residAuxStorage_.data_[3];

    const scalar_type dt_half = dt / static_cast< scalar_type >(2);
    const scalar_type t_phalf = t + dt_half;
    const scalar_type dt6 = dt / static_cast< scalar_type >( 6 );
    const scalar_type dt3 = dt / static_cast< scalar_type >( 3 );
    constexpr auto one  = ::pressio::utils::constants::one<scalar_type>();

    // stage 1: ytmp = y + auxRhs0*dt_half;
    policy_(y, auxRhs0, sys_.get(), t);
    ::pressio::containers::ops::do_update(ytmp, y, one, auxRhs0, dt_half);

    // stage 2: ytmp = y + auxRhs1*dt_half;
    policy_(ytmp, auxRhs1, sys_.get(), t_phalf);
    ::pressio::containers::ops::do_update(ytmp, y, one, auxRhs1, dt_half);

    // stage 3: ytmp = y + auxRhs2*dt;
    policy_(ytmp, auxRhs2, sys_.get(), t_phalf);
    ::pressio::containers::ops::do_update(ytmp, y, one, auxRhs2, dt);

    // stage 4: y_n += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    policy_(ytmp, auxRhs3, sys_.get(), t + dt);
    ::pressio::containers::ops::do_update(y, one,
					auxRhs0, dt6,
					auxRhs1, dt3,
					auxRhs2, dt3,
					auxRhs3, dt6);
  }//end doStep



  //  user-defined ops, with containers wrappers
  template<
    typename step_t,
    typename T = ops_t,
    typename _state_type = state_type,
    mpl::enable_if_t<
      std::is_void<T>::value == false and
      containers::meta::is_wrapper<_state_type>::value
      > * = nullptr
    >
  void doStep(_state_type & y,
  	      scalar_type t,
  	      scalar_type dt,
  	      step_t step){
    using op = typename ops_t::update_op;

    auto & ytmp	   = stateAuxStorage_.data_[0];
    auto & auxRhs0 = residAuxStorage_.data_[0];
    auto & auxRhs1 = residAuxStorage_.data_[1];
    auto & auxRhs2 = residAuxStorage_.data_[2];
    auto & auxRhs3 = residAuxStorage_.data_[3];

    const scalar_type dt_half = dt / static_cast< scalar_type >(2);
    const scalar_type t_phalf = t + dt_half;
    const scalar_type dt6 = dt / static_cast< scalar_type >( 6 );
    const scalar_type dt3 = dt / static_cast< scalar_type >( 3 );
    constexpr auto one  = ::pressio::utils::constants::one<scalar_type>();

    // stage 1: ytmp = y + auxRhs0*dt_half;
    policy_(y, auxRhs0, sys_.get(), t);
    op::do_update(*ytmp.data(), *y.data(), one, *auxRhs0.data(), dt_half);

    // stage 2: ytmp = y + auxRhs1*dt_half;
    policy_(ytmp, auxRhs1, sys_.get(), t_phalf);
    op::do_update(*ytmp.data(), *y.data(), one, *auxRhs1.data(), dt_half);

    // stage 3: ytmp = y + auxRhs2*dt;
    policy_(ytmp, auxRhs2, sys_.get(), t_phalf);
    op::do_update(*ytmp.data(), *y.data(), one, *auxRhs2.data(), dt);

    // stage 4: y_n += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    policy_(ytmp, auxRhs3, sys_.get(), t + dt);
    op::do_update(*y.data(), one,
    		  *auxRhs0.data(), dt6,
    		  *auxRhs1.data(), dt3,
    		  *auxRhs2.data(), dt3,
    		  *auxRhs3.data(), dt6);
  }//end doStep


#ifdef HAVE_PYBIND11
  /*
   * user does provide custom ops
   * interface with python
   */
  template<
    typename step_t,
    typename T = ops_t,
    typename _state_type = state_type,
    mpl::enable_if_t<
      std::is_void<T>::value == false and
      containers::meta::is_cstyle_array_pybind11<_state_type>::value
      > * = nullptr
    >
  void doStep(_state_type & y,
	      scalar_type t,
	      scalar_type dt,
	      step_t step){
    using op = typename ops_t::update_op;

    auto & ytmp	   = stateAuxStorage_.data_[0];
    auto & auxRhs0 = residAuxStorage_.data_[0];
    auto & auxRhs1 = residAuxStorage_.data_[1];
    auto & auxRhs2 = residAuxStorage_.data_[2];
    auto & auxRhs3 = residAuxStorage_.data_[3];

    const scalar_type dt_half = dt / static_cast< scalar_type >(2);
    const scalar_type t_phalf = t + dt_half;
    const scalar_type dt6 = dt / static_cast< scalar_type >( 6 );
    const scalar_type dt3 = dt / static_cast< scalar_type >( 3 );
    constexpr auto one  = ::pressio::utils::constants::one<scalar_type>();

    // stage 1: ytmp = y + auxRhs0*dt_half;
    policy_(y, auxRhs0, sys_.get(), t);
    op::do_update(ytmp, y, one, auxRhs0, dt_half);

    // stage 2: ytmp = y + auxRhs1*dt_half;
    policy_(ytmp, auxRhs1, sys_.get(), t_phalf);
    op::do_update(ytmp, y, one, auxRhs1, dt_half);

    // stage 3: ytmp = y + auxRhs2*dt;
    policy_(ytmp, auxRhs2, sys_.get(), t_phalf);
    op::do_update(ytmp, y, one, auxRhs2, dt);

    // stage 4: y_n += dt/6 * ( k1 + 2 * k2 + 2 * k3 + k4 )
    policy_(ytmp, auxRhs3, sys_.get(), t + dt);
    op::do_update(y, one,
    		  auxRhs0, dt6,
    		  auxRhs1, dt3,
    		  auxRhs2, dt3,
    		  auxRhs3, dt6);
  }//end doStep
#endif

}; //end class

}}}//end namespace pressio::ode::impl
#endif
