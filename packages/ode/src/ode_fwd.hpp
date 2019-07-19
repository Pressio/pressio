
#ifndef ODE_FORWARD_DECLARATIONS_HPP_
#define ODE_FORWARD_DECLARATIONS_HPP_

#include "ode_ConfigDefs.hpp"

namespace pressio{ namespace ode{

template<
  ExplicitEnum name,
  typename state_type,
  typename ...Args
  >
class ExplicitStepper;

template<
  ImplicitEnum name,
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename model_type,
  typename ...Args>
class ImplicitStepper;
//-----------------------------------

namespace policy{

template<
  typename state_type,
  typename model_type,
  typename velocity_type = state_type,
  typename enable = void>
class ExplicitVelocityStandardPolicy;

template<
  typename state_type,
  typename model_type,
  typename residual_type = state_type,
  typename enable = void>
class ImplicitResidualStandardPolicy;

template<
  typename state_type,
  typename model_type,
  typename jacobian_type,
  typename enable = void>
class ImplicitJacobianStandardPolicy;

#ifdef HAVE_PYBIND11
template<
  typename state_type,
  typename model_type,
  typename residual_type = state_type,
  typename enable = void>
class ImplicitResidualStandardPolicyPybind11;

template<
  typename state_type,
  typename model_type,
  typename jacobian_type,
  typename enable = void>
class ImplicitJacobianStandardPolicyPybind11;
#endif

}//end namespace policy
//-----------------------------------

namespace impl {

template<
  typename model_type,
  typename enable = void
  >
struct OdeSystemWrapper;

template<
  typename scalar_type,
  typename state_type,
  typename model_type,
  typename velocity_type,
  typename policy_type,
  typename ops,
  typename enable = void
  >
class ExplicitEulerStepperImpl;

template<
  typename scalar_type,
  typename state_type,
  typename model_type,
  typename velocity_type,
  typename policy_type,
  typename ops,
  typename enable = void
  >
class ExplicitRungeKutta4StepperImpl;

}//end namespace impl

}} // end namespace pressio::ode
#endif
