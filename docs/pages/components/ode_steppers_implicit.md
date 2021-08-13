
# ode: implicit steppers

Defined in: `<pressio/ode_steppers_implicit.hpp>`

Public namespace: `pressio::ode`


## Overview

pressio ode targets systems of the form:
@f[
\frac{d \boldsymbol{y}}{dt} =
\boldsymbol{f}(\boldsymbol{y},t; ...)
@f]

where @f$y@f$ is the state, @f$f@f$ is the RHS (also called velocity below), @f$t@f$ is time.
The ode package provides classes for doing time integration.
The design is rooted on the following main ideas:
- we separate the concept of a "stepper" from the actual "advancer": a stepper class
is responsible for knowing how to perform a step in time, while the "advancer" (or integrator)
is responbile of knowing how many times the stepper should be called.
- we use a policy-based design for evaluating the required operators,
such that this allows great flexiblity on using this time integration package.

The ode package is split into three main subcomponents:
- explicit methods
- implicit methods
- integrators

## Implicit Methods

### API

```cpp
template<class SystemType, class StateType>
auto create_keyword_stepper(const SystemType & system,
							const StateType & state);

template<class StateType, class ResidualPolicyType, class JacobianPolicyType>
auto create_keyword_stepper(const StateType & state,
							ResidualPolicyType && rPol,
							JacobianPolicyType && jPol);
```

where `keyword` is one of: `bdf1`, `bdf2`, `cranknicolson`.
