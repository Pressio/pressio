
# ode: explicit steppers

Defined in: `<pressio/ode_explicit.hpp>`

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


## Explicit Methods

### API:

```cpp
template<class SystemType, class StateType>
ReturnType create_keyword_stepper(const SystemType & system,
								  const StateType & state);
```

where `keyword` is one of: `forward_euler`, `runge_kutta4`, `adams_bashforth2`, `ssp_runge_kutta3`.


### Parameters

- `SystemType`:
  - class defining how to create an instance of the velocity @f$f@f$ and how to compute it;

- `StateType`:
  - type of the data structure you use for the state
  - Requirements: must be copy constructible


### Requirements

- The system class must conform to the following API:
  ```cpp
  struct ValidSystemForExplicitOde
  {
	using scalar_type   = /* */;
	using state_type    = /* */;
	using velocity_type = /* */;

	velocity_type createVelocity() const;
	void velocity(const state_type &, scalar_type time, velocity_type &) const;
  };
  ```

- the nested aliases `scalar_type`, `state_type` and `velocity_type` must be valid types:
these typedefs are needed for pressio to collect information

- if `StateType` is the type deduced for `state` from `create_...`, the following must hold:<br/>
  `std::is_same<StateType, typename ValidSystemForExplicitOde::state_type>::value == true`

