
# ode

Explicit methods are defined in: `<pressio/ode_explicit.hpp>`

Implicit methods are defined in: `<pressio/ode_implicit.hpp>`

Both can be includes via `<pressio/ode.hpp>`.

Public namespace: `pressio::ode`


## Overview
The ode package provides classes for doing time integration.
Specifically, the design of the package is rooted on the following main ideas:
- we separate the concept of a "stepper" from the actual "advancer": a stepper class
is responsible for knowing how to perform a step in time, while the "advancer" (or integrator)
is responbile of knowing how many times the stepper should be called.
- we use a policy-based design for evaluating the required operators,
such that this allows great flexiblity on using this time integration package.


## Explicit Methods

The explicit methods of the pressio ode aim to solve systems of the form:
@f[
\frac{d \boldsymbol{y}}{dt} =
\boldsymbol{f}(\boldsymbol{y},t; ...)
@f]

where @f$y@f$ is the state, @f$f@f$ is the RHS (also called velocity below), @f$t@f$ is time.

### API:
```cpp
template<class SystemType, class StateType>
ReturnType create_keyword_stepper(const SystemType & system,
								  const StateType & state);

template<class SystemType, class StateType, class policy_type>
ReturnType create_keyword_stepper(const SystemType & system,
							      const StateType & state,
                                  policy_type && rhsPolicy);
```
where `keyword` is one of: `forward_euler`, `runge_kutta4`, `adams_bashforth2`, `ssp_runge_kutta3`.


### Parameters and Requirements

- `SystemType`:
  - class defining the problem: how to create an instance of the velocity @f$f@f$ and how to compute it;
  - Requirements: must conform to the following API needed for doing explicit time integration
  ```cpp
  struct ValidSystemForExplicitOde
  {
	using scalar_type   = /* */;
	using state_type    = /* */;
	using velocity_type = /* typically the same as state_type */;

	velocity_type createVelocity() const;
	void velocity(const state_type &, scalar_type time, velocity_type &) const;
  };
  ```

- `StateType`:
  - type of the data structure you use for the state
  - Requirements: must be copy constructible

- `PolicyType`:
  - type of the class defining the policy for how to evaluate the right-hand-side
  - Requirements: must conform to a specific API
  ```cpp
  ```




## Implicit Methods
