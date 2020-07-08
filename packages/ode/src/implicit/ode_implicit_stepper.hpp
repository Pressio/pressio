
#ifndef ODE_IMPLICIT_METHODS_IMPLICIT_STEPPER_HPP_
#define ODE_IMPLICIT_METHODS_IMPLICIT_STEPPER_HPP_

#include "./impl/ode_implicit_stepper_composer_impl.hpp"

namespace pressio{ namespace ode{

template<typename stepper_tag, typename ...Args>
using ImplicitStepper =
  typename ::pressio::ode::implicitmethods::impl::Composer<stepper_tag, void, Args...>::type;

}} // end namespace pressio::ode
#endif
