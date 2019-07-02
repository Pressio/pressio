
#ifndef ODE_ENUM_HPP_
#define ODE_ENUM_HPP_

namespace pressio{ namespace ode{

enum class ExplicitEnum
  {Undefined, Euler, RungeKutta4};

enum class ImplicitEnum
  {Undefined, Euler, BDF2};

}}//end namespace pressio::ode
#endif
