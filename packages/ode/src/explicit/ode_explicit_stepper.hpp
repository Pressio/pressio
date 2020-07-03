
#ifndef ODE_EXPLICIT_METHODS_EXPLICIT_STEPPER_HPP_
#define ODE_EXPLICIT_METHODS_EXPLICIT_STEPPER_HPP_

#include "./impl/ode_explicit_stepper_impl.hpp"

namespace pressio{ namespace ode{ 

template<typename stepper_tag, typename ...Args>
using ExplicitStepper = typename ::pressio::ode::explicitmethods::impl::Stepper<stepper_tag, Args...>::type;

}} // end namespace pressio::ode
#endif
