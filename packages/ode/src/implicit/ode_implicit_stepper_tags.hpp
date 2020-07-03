
#ifndef ODE_IMPLICIT_STEPPER_TAGS_HPP_
#define ODE_IMPLICIT_STEPPER_TAGS_HPP_

namespace pressio{ namespace ode{ namespace implicitmethods{

struct Undefined{};
struct BDF1{};
struct BDF2{};

using Euler = BDF1;

// this is used to define a stepper that is user-defined,
// so there is no specific info about it and it can be
// anything since the user assembles the time-discrete operators
struct Arbitrary{};

}}} 
#endif
