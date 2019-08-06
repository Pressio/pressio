
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_EULER_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_EULER_STEPPER_IMPL_HPP_

#include "../ode_explicit_stepper_base.hpp"

namespace pressio{ namespace ode{ namespace impl{

template<
  typename scalar_type,
  typename state_type,
  typename system_type,
  typename velocity_type,
  typename residual_policy_type,
  typename ops_t
  >
class ExplicitEulerStepperImpl<scalar_type,
			       state_type,
			       system_type,
			       velocity_type,
			       residual_policy_type,
			       ops_t>
{

  static_assert( meta::is_legitimate_explicit_velocity_policy<
		 residual_policy_type>::value ||
		 meta::is_explicit_euler_velocity_standard_policy<
		 residual_policy_type>::value,
"EXPLICIT EULER VELOCITY_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");

  using this_t = ExplicitEulerStepperImpl< scalar_type,
					   state_type, system_type,
					   velocity_type,
					   residual_policy_type,
					   ops_t>;
  using resid_storage_t = OdeStorage<velocity_type, 1>;
  using system_wrapper_t = OdeSystemWrapper<system_type>;

  resid_storage_t residAuxStorage_;
  system_wrapper_t sys_;
  const residual_policy_type & policy_;

public:
  ExplicitEulerStepperImpl(const system_type & model,
			   const residual_policy_type & res_policy_obj,
			   const state_type & y0,
			   const velocity_type & r0)
    : residAuxStorage_(r0),
      sys_(model),
      policy_(res_policy_obj){}

  ExplicitEulerStepperImpl() = delete;
  ~ExplicitEulerStepperImpl() = default;

public:

  /*
   * user does NOT provide custom ops, so we use containers::ops
   */
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
    auto & auxRhs0 = residAuxStorage_.data_[0];
    //eval RHS
    policy_(y, auxRhs0, sys_.get(), t);
    // y = y + dt * rhs
    constexpr auto one  = ::pressio::utils::constants::one<scalar_type>();
    ::pressio::containers::ops::do_update(y, one, auxRhs0, dt);
  }

  /*
   * user does provide custom ops
   */
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
    auto & auxRhs0 = residAuxStorage_.data_[0];

    //eval RHS
    policy_(y, auxRhs0, sys_.get(), t);
    // y = y + dt * rhs
    constexpr auto one  = ::pressio::utils::constants::one<scalar_type>();
    op::do_update(*y.data(), one, *auxRhs0.data(), dt);
  }

#ifdef HAVE_PYBIND11
  /*
   * user does provide custom ops
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
    auto & auxRhs0 = residAuxStorage_.data_[0];

    policy_(y, auxRhs0, sys_.get(), t);
    // y = y + dt * rhs
    constexpr auto one  = ::pressio::utils::constants::one<scalar_type>();
    op::do_update(y, one, auxRhs0, dt);
  }
#endif

};

}}}//end namespace pressio::ode::impl
#endif
