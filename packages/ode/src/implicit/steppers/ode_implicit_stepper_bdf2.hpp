
#ifndef ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_BDF2_HPP_
#define ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_BDF2_HPP_

#include "ode_implicit_stepper_traits.hpp"
#include "ode_implicit_stepper_base.hpp"

namespace pressio{ namespace ode{

template<
  typename ode_state_type,
  typename ode_residual_type,
  typename ode_jacobian_type,
  typename system_type,
  typename ... Args
  >
class ImplicitStepper<
  ImplicitEnum::BDF2,
  ode_state_type,
  ode_residual_type,
  ode_jacobian_type,
  system_type,
  Args...
  >
  : public ImplicitStepperBase<
  ImplicitStepper<
    ImplicitEnum::BDF2,
    ode_state_type,
    ode_residual_type,
    ode_jacobian_type,
    system_type, Args...>,
  2 //num of aux states needed
  >
{
  using this_t	       = ImplicitStepper<ImplicitEnum::BDF2,
					 ode_state_type,
					 ode_residual_type,
					 ode_jacobian_type,
					 system_type,
					 Args...>;
  using stepper_base_t = ImplicitStepperBase<this_t, 2>;
  friend stepper_base_t;

  using mytraits       = details::traits<this_t>;
  using standard_res_policy_t = typename mytraits::standard_res_policy_t;
  using standard_jac_policy_t = typename mytraits::standard_jac_policy_t;
  using residual_pol_t = typename mytraits::residual_policy_t;
  using jacobian_pol_t = typename mytraits::jacobian_policy_t;
  using aux_stepper_t  = typename mytraits::aux_stepper_t;
  using scalar_t       = typename mytraits::scalar_t;
  static constexpr auto my_enum = mytraits::enum_id;

  aux_stepper_t & auxStepper_;

public:
  // these need to be here because are detected by solver
  using scalar_type	= scalar_t;
  using state_type	= ode_state_type;
  using residual_type	= ode_residual_type;
  using jacobian_type	= ode_jacobian_type;

public:
  ImplicitStepper() = delete;
  ~ImplicitStepper() = default;

  ImplicitStepper(const ode_state_type & y0,
  		  const system_type & model,
  		  const residual_pol_t & resPolicyObj,
  		  const jacobian_pol_t & jacPolicyObj,
		  aux_stepper_t & auxStepper)
    : stepper_base_t{y0, model, resPolicyObj, jacPolicyObj},
      auxStepper_{auxStepper}{}

  // cstr for standard residual and jacob policies
  template <
    typename T1 = standard_res_policy_t,
    typename T2 = standard_jac_policy_t,
    ::pressio::mpl::enable_if_t<
      mpl::is_same<T1, residual_pol_t>::value and
      mpl::is_same<T2, jacobian_pol_t>::value
      > * = nullptr
    >
  ImplicitStepper(const ode_state_type & y0,
		  const system_type & model,
		  aux_stepper_t & auxStepper)
    : stepper_base_t{y0, model},
      auxStepper_{auxStepper}{}

  // cstr for standard jacob policies
  template <
    typename T1 = standard_jac_policy_t,
    ::pressio::mpl::enable_if_t<
      mpl::is_same<T1, jacobian_pol_t>::value
      > * = nullptr
    >
  ImplicitStepper(const ode_state_type & y0,
  		  const system_type & model,
  		  const residual_pol_t & resPolicyObj,
		  aux_stepper_t & auxStepper)
    : stepper_base_t{y0, model, resPolicyObj},
      auxStepper_{auxStepper}{}

public:
  template<typename step_t,
	   typename solver_type>
  void operator()(ode_state_type & y,
		  scalar_t t,
		  scalar_t dt,
		  step_t step,
		  solver_type & solver){

    auto & auxY0 = this->stateAuxStorage_.data_[0];
    auto & auxY1 = this->stateAuxStorage_.data_[1];

    this->dt_ = dt;
    this->t_ = t;

    // first step, use auxiliary stepper
    if (step == 1){
      ::pressio::containers::ops::deep_copy(y, auxY0);
      auxStepper_(y, t, dt, step, solver);
    }
    if (step == 2){
      ::pressio::containers::ops::deep_copy(y, auxY1);
      solver.solve(*this, y);
    }
    if (step >= 3){
      ::pressio::containers::ops::deep_copy(auxY1, auxY0);
      ::pressio::containers::ops::deep_copy(y, auxY1);
      solver.solve(*this, y);
    }
  }

 template<typename step_t,
	  typename solver_type,
	  typename guess_callback_t>
  void operator()(ode_state_type & y,
		  scalar_t t,
		  scalar_t dt,
		  step_t step,
		  solver_type & solver,
		  guess_callback_t && guesserCb){

   auto & auxY0 = this->stateAuxStorage_.data_[0];
   auto & auxY1 = this->stateAuxStorage_.data_[1];

   this->dt_ = dt;
   this->t_ = t;

   // first step, use auxiliary stepper
   if (step == 1){
     ::pressio::containers::ops::deep_copy(y, auxY0);
     auxStepper_(y, t, dt, step, solver);
   }
   if (step == 2){
     ::pressio::containers::ops::deep_copy(y, auxY1);
     guesserCb(step, t, y);
     solver.solve(*this, y);
   }
   if (step >= 3){
     ::pressio::containers::ops::deep_copy(auxY1, auxY0);
     ::pressio::containers::ops::deep_copy(y, auxY1);
     guesserCb(step, t, y);
     solver.solve(*this, y);
   }
 }

};//end class

}} // end namespace pressio::ode
#endif
