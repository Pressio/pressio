
# ode: explicit steppers

Defined in: `<pressio/ode_steppers_explicit.hpp>`

Public namespace: `pressio::ode`


## Overview

Applicable to systems of the form:
@f[
\frac{d \boldsymbol{y}}{dt} =
\boldsymbol{f}(\boldsymbol{y},t; ...)
@f]

where @f$y@f$ is the state, @f$f@f$ is the RHS (also called velocity below), @f$t@f$ is time.<br/>
Explicit methods calculate the state of a system at a later time
from the state of the system at the current time and potentially previous times.

## API:

```cpp
template<class StateType, class SystemType>
auto create_keyword_stepper(const StateType & state,
	                        const SystemType & system);
```

where `keyword` is one of: `forward_euler`, `runge_kutta4`, `adams_bashforth2`, `ssp_runge_kutta3`.

## Parameters

- `StateType`:
  - type of the data structure you use for the state

- `SystemType`:
  - class defining how to create an instance of the velocity @f$f@f$ and how to compute it;


## Requirements

- `StateType`: must be copy constructible

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
these typedefs are detected by pressio

- if `StateType` is the type deduced for `state` from `create_...`, the following must hold:<br/>
  `std::is_same<StateType, typename ValidSystemForExplicitOde::state_type>::value == true`

- the create function above returns an instance of the desired stepper.
  The stepper object returned satisfies the "steppable" concept discussed [here](/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio/docs/html/md_pages_components_ode_advance.html), so one can use the "advancers" functions to step forward.

## Ops

When using custom data types not supported in [pressio ops](/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio/docs/html/md_pages_components_ops.html), you need to specialize a trait class and some operations
such that pressio can operate on your data. For the sake of explanation, suppose that you use:

```cpp
using scalar_type   = double;
using state_type    = ACustomStateType;  //this can be any type
```

Then you need to provide the following specializations:

```cpp
namespace pressio{

template<> struct Traits<ACustomStateType>{
  using scalar_type = double;
};

namespace ops{

ACustomStateType clone(const ACustomStateType & src){
  /* return a deep copy of src */
}

void set_zero(ACustomStateType & object){
  /* set elements zero */
}

void update(ACustomStateType & v,
            const ACustomStateType & v1, const scalar_type b)
{
  // elementwise compute : v = b*v1
}

void update(ACustomStateType & v,        const scalar_type a,
		    const ACustomStateType & v1, const scalar_type b)
{
  // elementwise compute : v = a*v + b*v1
}

void update(ACustomStateType & v,
            const ACustomStateType & v0, const scalar_type a,
            const ACustomStateType & v1, const scalar_type b)
{
  // elementwise compute: v = a*v0 + b*v1
}

void update(ACustomStateType & v,        const scalar_type c,
            const ACustomStateType & v0, const scalar_type a,
            const ACustomStateType & v1, const scalar_type b)
{
  // elementwise compute : v = c*v + a*v0 + b*v1
}

void update(ACustomStateType & v,
			const ACustomStateType & v1, const scalar_type b,
			const ACustomStateType & v2, const scalar_type c,
			const ACustomStateType & v3, const scalar_type d)
{
  // elementwise compute: v = b*v1 + c*v2 + d*v3
}

void update(ACustomStateType & v,		 const scalar_type a,
			const ACustomStateType & v1, const scalar_type b,
			const ACustomStateType & v2, const scalar_type c,
			const ACustomStateType & v3, const scalar_type d)
{
  // elementwise compute: v = a*v + b*v1 + c*v2 + d*v3
}

void update(ACustomStateType & v,
			const ACustomStateType & v1, const scalar_type b,
			const ACustomStateType & v2, const scalar_type c,
			const ACustomStateType & v3, const scalar_type d,
			const ACustomStateType & v4, const scalar_type e)
{
  // elementwise compute: v = b*v1 + c*v2 + d*v3 + e*v4
}

void update(ACustomStateType & v,		 const scalar_type a,
			const ACustomStateType & v1, const scalar_type b,
			const ACustomStateType & v2, const scalar_type c,
			const ACustomStateType & v3, const scalar_type d,
			const ACustomStateType & v4, const scalar_type e)
{
  // elementwise compute: v = a*v + b*v1 + c*v2 + d*v3 + e*v4
}

}}//end namepsace pressio::ops
```

## Example usage
```cpp
int main()
{
  // assuming that:
  // stateObj  is the state
  // systemObj is the system instance

  namespace pode = pressio::ode;
  auto stepper = pode::create_forward_euler_stepper(stateObj, systemObj);

  // use the stepper to advance in time,
  // for example using the advancer function
  const double time0 = 0.;
  const double dt = 0.1;
  const pode::step_count_type num_steps = 100;
  pode::advance_n_steps(stepper, stateObj, time0, dt, num_steps);
}
```
