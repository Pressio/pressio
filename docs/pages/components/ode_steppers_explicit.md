
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
In pressio, a "stepper" is an abstraction that represents the "how" to take a step.


## API

```cpp
template<class StateType, class SystemType>
auto create_keyword_stepper(const StateType & state,
	                        const SystemType & system);
```

where `keyword` is one of: `forward_euler`, `runge_kutta4`, `adams_bashforth2`, `ssp_runge_kutta3`.
This function returns an instance of the desired stepper.
The returned stepper object satisfies the "steppable" concept discussed [here](/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio/docs/html/md_pages_components_ode_advance.html), so one can use the "advancers" functions to step forward.


## Parameters

- `StateType`:
  - type of the data structure you use for the state

- `SystemType`:
  - class defining how to create an instance of the velocity @f$f@f$ and how to compute it;


## Requirements

- `StateType`: must be copy constructible

- The system class must conform to the following API:

  ```cpp
  struct SystemForExplicitOde
  {
	using scalar_type   = /* */;
	using state_type    = /* */;
	using velocity_type = /* */;

	velocity_type createVelocity() const;
	void velocity(const state_type &, scalar_type time, velocity_type &) const;
  };
  ```

  the nested aliases `scalar_type`, `state_type` and `velocity_type` must be *valid* types since
  they are detected by pressio

- if `StateType` is the type deduced for `state` from `create_...`, the following must hold:<br/>
  `std::is_same<StateType, typename SystemForExplicitOde::state_type>::value == true`

## Example usage

```cpp
#include "pressio/type_traits.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressio/ode_steppers_explicit.hpp"
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


## Required specializations for custom types

When using custom data types not supported in [pressio ops](/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio/docs/html/md_pages_components_ops.html), you need to provide specializations of a trait class and certain operations
and make them "visible" to the compiler to find them and such that pressio can operate on your data.
For the sake of explanation, suppose that you use `double`
as value type and `ACustomStateType` is what you use for the state, then you would need to do something like this:

```cpp
#include "pressio/type_traits.hpp"

// assuming ACustomStateType has already been declared

namespace pressio{

template<> struct Traits<ACustomStateType>{
  using scalar_type = double;
};

namespace ops{

void deep_copy(ACustomStateType & dest, const ACustomStateType & src){
  /* deep copy src into dest */
}

ACustomStateType clone(const ACustomStateType & src){
  /* return a deep copy of src */
}

void set_zero(ACustomStateType & object){
  /* set elements to zero */
}

void update(ACustomStateType & v,        const double a,
		    const ACustomStateType & v1, const double b)
{
  // elementwise compute : v = a*v + b*v1
}

void update(ACustomStateType & v,        const double a,
            const ACustomStateType & v1, const double b,
            const ACustomStateType & v2, const double c)
{
  // elementwise compute : v = a*v + b*v1 + c*v2
}

void update(ACustomStateType & v,		 const double a,
			const ACustomStateType & v1, const double b,
			const ACustomStateType & v2, const double c,
			const ACustomStateType & v3, const double d)
{
  // elementwise compute: v = a*v + b*v1 + c*v2 + d*v3
}

void update(ACustomStateType & v,		 const double a,
			const ACustomStateType & v1, const double b,
			const ACustomStateType & v2, const double c,
			const ACustomStateType & v3, const double d,
			const ACustomStateType & v4, const double e)
{
  // elementwise compute: v = a*v + b*v1 + c*v2 + d*v3 + e*v4
}
}}//end namepsace pressio::ops

#include "pressio/ode_advancers.hpp"
#include "pressio/ode_steppers_explicit.hpp"

int main()
{
  // same code as shown above
}
```

Note that in the snippet above the order of the include statements matter:
this is because your `Trait` and kernel specializations need to be found by the compiler.
However, to make the code cleaner, you can obviously move all kernels specializations
to a separate header file, but make sure to keep the correct order, for example as follows:

```cpp
#include "pressio/type_traits.hpp"
#include "my_specializations.hpp" // contains all your specializations
#include "pressio/ode_advancers.hpp"
#include "pressio/ode_steppers_explicit.hpp"
int main()
{
  // same code as shown above
}
```
