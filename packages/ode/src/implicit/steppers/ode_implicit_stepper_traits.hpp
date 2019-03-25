
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_TRAITS_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_TRAITS_HPP_

#include "../../ode_forward_declarations.hpp"

namespace rompp{ namespace ode{

namespace impl{

template <ImplicitEnum whichone>
struct orderHelper;

template <>
struct orderHelper<ImplicitEnum::Euler>{
  using order_t	= unsigned int;
  static constexpr order_t order_value = 1;
};

template <>
struct orderHelper<ImplicitEnum::BDF2>{
  using order_t	= unsigned int;
  static constexpr order_t order_value = 2;
};
//--------------------------------------------

template <ImplicitEnum whichone>
struct stepsHelper;

template <>
struct stepsHelper<ImplicitEnum::Euler>{
  static constexpr unsigned int steps = 1;
};

template <>
struct stepsHelper<ImplicitEnum::BDF2>{
  static constexpr unsigned int steps = 2;
};
//--------------------------------------------



template <
  typename state_type, typename residual_type,
  typename jacobian_type, typename model_type
  >
struct standardPolHelper{
  using residual_policy_t =
    policy::ImplicitResidualStandardPolicy<state_type, model_type, residual_type>;
  using jacobian_policy_t =
    policy::ImplicitJacobianStandardPolicy<state_type, model_type, jacobian_type>;
};

template <
  typename state_type, typename jacobian_type, typename model_type
  >
struct standardPolHelper<state_type, void, jacobian_type, model_type>
{
  using jacobian_policy_t =
    policy::ImplicitJacobianStandardPolicy<state_type, model_type, jacobian_type>;
};

template <
  typename state_type, typename residual_type, typename model_type
  >
struct standardPolHelper<state_type, residual_type, void, model_type>
{
  using residual_policy_t =
    policy::ImplicitResidualStandardPolicy<state_type, model_type, residual_type>;
};

//--------------------------------------------




// metaf to help define the implementation class type
template <
  ImplicitEnum whicone,
  typename state_type, typename residual_type,
  typename jacobian_type, typename model,
  typename aux_stepper_type,  typename res_pol_t,
  typename jac_pol_t, typename enable = void
  >
struct implClassHelper;



//******************
//****** EULER *****
//******************

// EULER with STANDARD policies
template <
  typename state_type, typename residual_type,
  typename jacobian_type, typename model
  >
struct implClassHelper<
  ImplicitEnum::Euler, state_type, residual_type,
  jacobian_type, model, void, void, void
  > : standardPolHelper<state_type, residual_type, jacobian_type, model>
{
  using base_t = standardPolHelper<state_type, residual_type, jacobian_type, model>;
  using typename base_t::residual_policy_t;
  using typename base_t::jacobian_policy_t;
  using impl_t = impl::ImplicitEulerStepperImpl
    <state_type, residual_type, jacobian_type,
     model, residual_policy_t, jacobian_policy_t>;
};


// EULER with:
// residual = user-defined
// jacobian = standard
template <
  typename state_type, typename residual_type,
  typename jacobian_type, typename model,
  typename res_pol_t
  >
struct implClassHelper<
  ImplicitEnum::Euler, state_type, residual_type,
  jacobian_type, model, void, res_pol_t, void,
  core::meta::enable_if_t<
    !std::is_void<res_pol_t>::value
    >
  > : standardPolHelper<state_type, void, jacobian_type, model>
{
  using base_t = standardPolHelper<state_type, void, jacobian_type, model>;
  using typename base_t::jacobian_policy_t;

  using residual_policy_t = res_pol_t;

  using impl_t = impl::ImplicitEulerStepperImpl
    <state_type, residual_type, jacobian_type,
     model, residual_policy_t, jacobian_policy_t>;
};


// EULER and user-defined policies
template <
  typename state_type, typename residual_type,
  typename jacobian_type, typename model,
  typename res_pol_t, typename jac_pol_t
  >
struct implClassHelper<
  ImplicitEnum::Euler, state_type, residual_type,
  jacobian_type, model, void, res_pol_t, jac_pol_t,
  core::meta::enable_if_t<
    !std::is_void<res_pol_t>::value and !std::is_void<jac_pol_t>::value
    >
  >
{
  using residual_policy_t = res_pol_t;
  using jacobian_policy_t = jac_pol_t;
  using impl_t = impl::ImplicitEulerStepperImpl<state_type, residual_type, jacobian_type,
						model,res_pol_t, jac_pol_t>;
};



//******************
//****** BDF2 *****
//******************

// BDF2 and STANDARD policies
template <
  typename state_type, typename residual_type,
  typename jacobian_type, typename model,
  typename aux_step_t
  >
struct implClassHelper<
  ImplicitEnum::BDF2, state_type, residual_type,
  jacobian_type, model, aux_step_t, void, void
  > : standardPolHelper<state_type, residual_type, jacobian_type, model>
{
  using base_t = standardPolHelper<state_type, residual_type, jacobian_type, model>;
  using typename base_t::residual_policy_t;
  using typename base_t::jacobian_policy_t;
  using impl_t = impl::ImplicitBDF2StepperImpl
    <state_type, residual_type, jacobian_type, model,
     aux_step_t, residual_policy_t, jacobian_policy_t>;
};


// BDF2 with:
// residual = user-defined
// jacobian = standard
template <
  typename state_type, typename residual_type,
  typename jacobian_type, typename model,
  typename aux_step_t, typename res_pol_t
  >
struct implClassHelper<
  ImplicitEnum::BDF2, state_type, residual_type,
  jacobian_type, model, aux_step_t, res_pol_t, void,
  core::meta::enable_if_t<
    !std::is_void<res_pol_t>::value
    >
  > : standardPolHelper<state_type, void, jacobian_type, model>
{
  using base_t = standardPolHelper<state_type, void, jacobian_type, model>;
  using typename base_t::jacobian_policy_t;

  using residual_policy_t = res_pol_t;

  using impl_t = impl::ImplicitBDF2StepperImpl
    <state_type, residual_type, jacobian_type,
     model, aux_step_t, residual_policy_t, jacobian_policy_t>;
};


// BDF2 and user-defined policies
template <
  typename state_type, typename residual_type,
  typename jacobian_type, typename model, typename aux_step_t,
  typename res_pol_t, typename jac_pol_t
  >
struct implClassHelper<
  ImplicitEnum::BDF2, state_type, residual_type,
  jacobian_type, model, aux_step_t, res_pol_t, jac_pol_t,
  core::meta::enable_if_t<
    !std::is_void<res_pol_t>::value and
    !std::is_void<jac_pol_t>::value
    >
  >
{
  using residual_policy_t = res_pol_t;
  using jacobian_policy_t = jac_pol_t;
  using impl_t = impl::ImplicitBDF2StepperImpl
    <state_type, residual_type, jacobian_type,
     model, aux_step_t, res_pol_t, jac_pol_t>;
};

}//end namespace impl



namespace details{

// traits class suitable for any implicit stepper
template<
  ImplicitEnum whichone,
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename model_type,
  typename aux_stepper_type,
  typename residual_policy_type,
  typename jacobian_policy_type
  >
struct traits<
  ImplicitStepper<
    whichone,
    state_type, residual_type, jacobian_type,
    model_type, aux_stepper_type,
    residual_policy_type, jacobian_policy_type
    >
  > : impl::orderHelper<whichone>,
      impl::stepsHelper<whichone>,
      impl::implClassHelper<whichone, state_type, residual_type, jacobian_type,
			    model_type, aux_stepper_type,
			    residual_policy_type, jacobian_policy_type>
{

  static constexpr ode::ImplicitEnum enum_id = whichone;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;

  using state_t		  = state_type;
  using residual_t	  = residual_type;
  using jacobian_t	  = jacobian_type;
  using model_t		  = model_type;
  using aux_stepper_t	  = aux_stepper_type;
  using scalar_t	  = typename core::details::traits<state_type>::scalar_t;

  using impl::orderHelper<whichone>::order_value;
  using impl::stepsHelper<whichone>::steps;
  using impl_help_t = impl::implClassHelper<whichone, state_type, residual_type,
					    jacobian_t, model_t, aux_stepper_t,
					    residual_policy_type, jacobian_policy_type>;
  using typename impl_help_t::impl_t;
  using typename impl_help_t::residual_policy_t;
  using typename impl_help_t::jacobian_policy_t;

};

}}}//end namespace rompp::ode::details
#endif
