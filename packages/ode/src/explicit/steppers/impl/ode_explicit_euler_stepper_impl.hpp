
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_EULER_STEPPER_IMPL_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_IMPL_EXPLICIT_EULER_STEPPER_IMPL_HPP_

#include "../ode_explicit_stepper_base.hpp"

namespace rompp{ namespace ode{ namespace impl{

template<
  typename scalar_type,
  typename ode_state_type,
  typename model_type,
  typename ode_residual_type,
  typename residual_policy_type,
  typename ops_t
  >
class ExplicitEulerStepperImpl<scalar_type,
			       ode_state_type,
			       model_type,
			       ode_residual_type,
			       residual_policy_type,
			       ops_t>
  : private OdeStorage<ode_state_type, ode_residual_type, 0, 1>,
    private ExplicitOdeAuxData<model_type, residual_policy_type>
{

  static_assert( meta::is_legitimate_explicit_residual_policy<
		 residual_policy_type>::value ||
		 meta::is_explicit_euler_residual_standard_policy<
		 residual_policy_type>::value,
"EXPLICIT EULER RESIDUAL_POLICY NOT ADMISSIBLE, \
MAYBE NOT A CHILD OF ITS BASE OR DERIVING FROM WRONG BASE");

  using this_t = ExplicitEulerStepperImpl< scalar_type,
					   ode_state_type, model_type,
					   ode_residual_type,
					   residual_policy_type,
					   ops_t>;
  using storage_base_t = OdeStorage<ode_state_type, ode_residual_type, 0, 1>;
  using auxdata_base_t = ExplicitOdeAuxData<model_type, residual_policy_type>;

  using storage_base_t::auxRHS_;
  using auxdata_base_t::model_;
  using auxdata_base_t::residual_obj_;

public:
  ExplicitEulerStepperImpl(const model_type & model,
			   const residual_policy_type & res_policy_obj,
			   const ode_state_type & y0,
			   const ode_residual_type & r0)
    : storage_base_t(r0), auxdata_base_t(model, res_policy_obj){}

  ExplicitEulerStepperImpl(const residual_policy_type & res_policy_obj,
			   const ode_state_type & y0,
			   const ode_residual_type & r0)
    : storage_base_t(r0), auxdata_base_t(res_policy_obj){}

  ExplicitEulerStepperImpl() = delete;
  ~ExplicitEulerStepperImpl() = default;

public:

  // if user does NOT provide ops, then user wrapper ops
  template<typename step_t,
	   typename T = ops_t,
	   mpl::enable_if_t<
	     std::is_void<T>::value
	     > * = nullptr
	   >
  void doStep(ode_state_type & y,
	      scalar_type t,
	      scalar_type dt,
	      step_t step){
    //eval RHS
    this->evalRHS(y, auxRHS_[0], t);
    // y = y + dt * rhs
    constexpr auto one  = ::rompp::core::constants::one<scalar_type>();
    ::rompp::core::ops::do_update(y, one, auxRHS_[0], dt);
  }

  // if user provides ops, then use them
  template<typename step_t,
	   typename T = ops_t,
	   mpl::enable_if_t<
	     std::is_void<T>::value == false
	     > * = nullptr
	   >
  void doStep(ode_state_type & y,
	      scalar_type t,
	      scalar_type dt,
	      step_t step){
    using op = typename ops_t::update_op;
    //eval RHS
    this->evalRHS(y, auxRHS_[0], t);
    // y = y + dt * rhs
    constexpr auto one  = ::rompp::core::constants::one<scalar_type>();
    op::do_update(*y.data(), one, *auxRHS_[0].data(), dt);
  }

private:
  template<typename T = model_type,
	   typename std::enable_if<
	     std::is_void<T>::value
	     >::type * = nullptr>
  void evalRHS(const ode_state_type & y,
	       ode_residual_type & RHS,
	       scalar_type t){
    (*residual_obj_)(y, RHS, t);
  }

  template<typename T = model_type,
	   typename std::enable_if<
	     !std::is_void<T>::value
	     >::type * = nullptr>
  void evalRHS(const ode_state_type & y,
	       ode_residual_type & RHS,
	       scalar_type t){
    (*residual_obj_)(y, RHS, *model_, t);
  }
};

}}}//end namespace rompp::ode::impl
#endif
