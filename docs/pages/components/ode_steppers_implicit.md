
# ode: implicit steppers

Defined in: `<pressio/ode_steppers_implicit.hpp>`

Public namespace: `pressio::ode`

## Overview

Represents the concept of a stepper in the context of implicit methods.
Recall that implicit methods calculate the state of a system at the next time
by solving a system of equations involving both the current state of
the system and the later one.

This package in pressio lets you create a stepper
for any system written in *continuous-time* form:
@f[
\frac{d \boldsymbol{y}}{dt} =
\boldsymbol{f}(\boldsymbol{y},t; ...)
@f]

and/or in a *discrete-time* form
@f[
\boldsymbol{R}(\boldsymbol{y}, \boldsymbol{y_{n-1}}, ..., t_n, dt_n; ...) = \boldsymbol{0}
@f]

Here, @f$y@f$ is the full-order model (FOM) state,
@f$f@f$ the FOM velocity, @f$t@f$ is time, and @f$R@f$ is the residual.


## API for continuous-time systems

```cpp
// overload 1
template<class StateType, class SystemType>
auto create_keyword_stepper(const StateType & state,
							const SystemType & system);

// overload 2
template<class StateType, class ResidualPolicyType, class JacobianPolicyType>
auto create_keyword_stepper(const StateType & state,
							ResidualPolicyType && rPol,
							JacobianPolicyType && jPol);
```

where:

| `keyword`         	| Method                  	| Discrete Residual Formula                                                                                      	|
|-----------------	|-------------------------	|----------------------------------------------------------------------------------------------	|
| `bdf1`          	| Backward Diff 1st order 	| @f$R = y_{n+1}-y_{n}- hf(t_{n+1},y_{n+1})@f$                                                      	|
| `bdf2`          	| Backward Diff 2nd order 	| @f$R = y_{n+1}-{\tfrac {4}{3}}y_{n}+{\tfrac {1}{3}}y_{n-1} - {\tfrac {2}{3}}hf(t_{n+1},y_{n+1})@f$ 	|
| `cranknicolson` 	| Crank-Nicolson          	| @f$R = y_{n+1}- y_{n} - {\tfrac {1}{2}} h \left( f(t_{n+1},y_{n+1}) + f(t_{n},y_{n}) \right)@f$  	|

## Parameters

- `StateType`:
  - type of the data structure you use for the state

- `SystemType`:
  - class defining how to create an instance of the velocity @f$f@f$ and how to compute it;

- `ResidualPolicyType`, `JacobianPolicyType`:
  - policy types if you want to use custom ones to compute the discrete operators


Notes:
- if you use the first overload above, pressio uses default policies
  to compute the residual and jacobian.

- the second overload allows you to define custom policies to compute
  the discrete operators


## Requirements

- `StateType`: must be copy constructible

- The system class must conform to the following API:
  ```cpp
  struct SystemForImplicitOde
  {
	using scalar_type   = /* */;
	using state_type    = /* */;
	using velocity_type = /* */;
	using jacobian_type =  /* */;

	velocity_type createVelocity() const;
    jacobian_type createJacobian() const;
	void velocity(const state_type &, scalar_type time, velocity_type &) const;
    void jacobian(const state_type &, scalar_time time, jacobian_type &) const;
  };
  ```

  the nested type aliases must be *valid* types since they are detected by pressio

- if `StateType` is the type deduced for `state` from `create_...`, the following must hold:<br/>
  `std::is_same<StateType, typename SystemForImplicitOde::state_type>::value == true`

- `ResidualPolicyType` must conform to:
  ```cpp
  class ResidualPolicy
  {
  public:
	using residual_type = /* this type alias needs to be found */;

	residual_type create() const;

	template <class TagType, class AuxStatesType, class AuxRhsType, class TimeType, class StepType>
	void compute(const StateType & y,
				 const AuxStatesType & auxStates,
				 AuxRhsType & auxRhs,
				 const TimeType & t,
				 const TimeType & dt,
				 StepType step,
				 residual_type & R) const;
  };
  ```

- `JacobianPolicyType`:
  *todo*


### If you use custom policies:

If you want to use custom policies for computing residual and Jacobian,
you need are responsible for ensuring things are correct.

| TagType         	| Auxiliary States contain: 	| Auxiliary RHS contain:  	|
|-----------------	|---------------------------	|------------------------	|
| `BDF1`          	| `y_n`                     	| (no velocities stored) 	|
| `BDF2`          	| `y_n`, `y_{n-1}`          	| (no velocities stored) 	|
| `CrankNicolson` 	| `y_n`                     	| `f_n`, `f_{n+1}`        	|



## API for discrete-time systems

```cpp
template<int num_states, class StateType, class SystemType>
auto create_arbitrary_stepper(const StateType & state,
						      SystemType && system);
```
