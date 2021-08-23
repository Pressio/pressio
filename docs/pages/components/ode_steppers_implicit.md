
# ode: implicit steppers

Defined in: `<pressio/ode_steppers_implicit.hpp>`

Public namespace: `pressio::ode`

## Overview

Provides functionalities to create a stepper suitable for implicit methods.
Recall that implicit methods calculate the state of a system at the next time
by solving a system of equations involving both the current state of
the system and the later one. An implicit stepper is an object that knows
how to take one *implicit* step.

Pressio implicit steppers are applicable to any system written in *continuous-time* form:
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

@code{.cpp}
// overload 1
template<class StateType, class SystemType>
auto create_keyword_stepper(const StateType & state,
							const SystemType & system);

// overload 2
template<class StateType, class ResidualPolicyType, class JacobianPolicyType>
auto create_keyword_stepper(const StateType & state,
							ResidualPolicyType && rPol,
							JacobianPolicyType && jPol);
@endcode

where the currently choices are:

| `keyword`         	| Method                  	| Discrete Residual Formula                                                                                      	|
|-----------------	|-------------------------	|----------------------------------------------------------------------------------------------	|
| `bdf1`          	| Backward Diff 1st order 	| @f$R = y_{n+1}-y_{n}- hf(t_{n+1},y_{n+1})@f$                                                      	|
| `bdf2`          	| Backward Diff 2nd order 	| @f$R = y_{n+1}-{\tfrac {4}{3}}y_{n}+{\tfrac {1}{3}}y_{n-1} - {\tfrac {2}{3}}hf(t_{n+1},y_{n+1})@f$ 	|
| `cranknicolson` 	| Crank-Nicolson          	| @f$R = y_{n+1}- y_{n} - {\tfrac {1}{2}} h \left( f(t_{n+1},y_{n+1}) + f(t_{n},y_{n}) \right)@f$  	|

The create functions above return an instance of the desired stepper,
which satisfies the "steppable" concept discussed [here](/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio/docs/html/md_pages_components_ode_advance.html). Therefore, one can use
the "advancers" functions to step forward.


### Parameters

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


### Requirements

- `StateType`: must be copy constructible

- `SystemType` must conform to the following API:
  @code{.cpp}
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
  @endcode

  the nested type aliases must be *valid* types since they are detected by pressio

- if `StateType` is the type deduced for `state` passed to `create_...`, the following must hold:<br/>
  `std::is_same<StateType, typename SystemForImplicitOde::state_type>::value == true`

- `ResidualPolicyType` must conform to:
  @code{.cpp}
  class ResidualPolicy
  {
  public:
	using residual_type = /* this type alias needs to be found */;

	residual_type create() const;

	template <class TagType, class AuxStatesType, class AuxRhsType, class TimeType, class StepType>
	void compute(const StateType & state,
				 const AuxStatesType & auxStates,
				 AuxRhsType & auxRhs,
				 const TimeType & time_at_n_plus_one,
				 const TimeType & dt,
				 StepType step,
				 residual_type & R) const;
  };
  @endcode

- `JacobianPolicyType`:
  @code{.cpp}
  class JacobianPolicy
  {
  public:
	using jacobian_type = /* this type alias needs to be found */;

	jacobian_type create() const;

	template <class TagType, class AuxStatesType, class TimeType, class StepType>
	void compute(const StateType & state,
				 const AuxStatesType & auxStates,
				 const TimeType & time_at_n_plus_one,
				 const TimeType & dt,
				 StepType step,
				 jacobian_type & J) const;
  };
  @endcode


<!-- ### Implicit stepper class synopsis -->
<!-- @code{.cpp} -->
<!-- class Stepper -->
<!-- { -->
<!-- public: -->
<!--   // these aliases are detected by solver -->
<!--   using scalar_type	= ScalarType; -->
<!--   using state_type	= StateType; -->
<!--   using residual_type	= ResidualType; -->
<!--   using jacobian_type	= JacobianType; -->

<!--   using tag_name = /* tag identifying the scheme */ -->
<!--   static constexpr bool is_implicit = true; -->
<!--   static constexpr bool is_explicit = false; -->

<!--   template<typename SolverType, typename ...Args> -->
<!--   void operator()(state_type & odeState, -->
<!-- 				  const ScalarType & currentTime, -->
<!-- 				  const ScalarType & dt, -->
<!-- 				  const int32_t & stepNumber, -->
<!-- 				  SolverType & solver, -->
<!-- 				  Args&& ...args); -->

<!--   ResidualType createResidual() const; -->
<!--   JacobianType createJacobian() const; -->
<!--   void residual(const StateType & odeState, ResidualType & R) const; -->
<!--   void jacobian(const StateType & odeState, JacobianType & J) const; -->
<!-- }; -->
<!-- @endcode -->


### If you use custom policies:

If you want to use custom policies for computing residual and Jacobian,
you need are responsible for ensuring things are correct.
In particular, you should be aware of the following:

- `State`: argument passed to `compute` method contains the prediction at `n+1`.

- `auxStates`, `auxRhs` (whose type you don't need to know) contain the needed auxiliary
states and RHS evaluations, respectively, needed to compute the operators
for a certain scheme. All you need to know about these containers is the following:


| TagType         	| Description/Info 	|
|-----------------	|---------------------------	|
| `BDF1`          	| `auxStates`: contains: state at n-th step <br/> &emsp; &emsp; &emsp; &emsp;  &ensp; Use: `const auto & yn = auxStates(pressio::ode::n());` <br/> `auxRhs`: Empty |
| `BDF2`          	| `auxStates`: contains: states at n-th and (n-1)-th step <br/> &emsp; &emsp; &emsp; &emsp;  &ensp; Use: `const auto & yn = auxStates(pressio::ode::n());` <br/> &emsp; &emsp; &emsp; &emsp; &ensp; `const auto & ynm1 = auxStates(pressio::ode::nMinusOne());` <br/> `auxRhs`: Empty |
| `CrankNicolson` | `auxStates`: contains: states at n-th step <br/>  &emsp; &emsp; &emsp; &emsp; &ensp; Use: `const auto & yn = auxStates(pressio::ode::n());` <br/> `auxRhs`: contains evaluations of the RHS are n-th and (n+1)-th steps <br/>  &emsp; &emsp; &emsp; Use: `auto & fn = auxRhs(pressio::ode::n());` <br/> &ensp; &ensp; &ensp; &emsp; `auto & fnp1 = auxRhs(pressio::ode::nPlusOne());` |


### Example usage

@code{.cpp}
#include "pressio/type_traits.hpp"
#include "pressio/ode_solvers_nonlinear.hpp"
#include "pressio/ode_advancers.hpp"
#include "pressio/ode_steppers_implicit.hpp"
int main()
{
  // assuming that:
  // stateObj  is the state
  // systemObj is the system instance

  namespace pode = pressio::ode;
  auto stepper = pode::create_bdf1_stepper(stateObj, systemObj);

  // create a solver, here for simplicity we show the case where
  // for the types used, we can leverage pressio solvers
  using jacobian_t = typename problem_t::jacobian_type;
  using lin_solver_t = pressio::linearsolvers::Solver</*some tag to specify method*/, jacobian_t>;
  lin_solver_t linSolverObj;
  auto nonLinSolver = nonlinearsolvers::create_newton_raphson(stepperObj,y,linSolverObj);

  // use the stepper to advance in time,
  // for example using the advancer function
  const double time0 = 0.;
  const double dt = 0.1;
  const pode::step_count_type num_steps = 100;
  pode::advance_n_steps(stepper, stateObj, time0, dt, num_steps, nonLinearSolver);
}
@endcode


### Required specializations for custom types

When using custom data types not supported in [pressio ops](/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio/docs/html/md_pages_components_ops.html), you need to provide specializations of a trait class and certain operations
and make them "visible" to the compiler to find them and such that pressio can operate on your data.
For the sake of explanation, suppose that you use `double`
as value type and `ACustomStateType` is what you use for the state, `ACustomMatrixType` is what you
use for matrix, then you would need to do something like this:

@code{.cpp}
#include "pressio/type_traits.hpp"

// assuming ACustomStateType has already been declared
// assuming ACustomMatrixType has already been declared

namespace pressio{

template<> struct Traits<ACustomStateType>{
  using scalar_type = double;
};

template<> struct Traits<ACustomMatrixType>{
  using scalar_type = double;
};

namespace ops{

void deep_copy(ACustomStateType & dest, const ACustomStateType & src){
  /* deep copy src into dest */
}

ACustomStateType clone(const ACustomStateType & src){
  /* return a deep copy of src */
}

void scale(ACustomMatrixType & M, double factor){
  /* scale elementwise by factor */
}

void add_to_diagonal(ACustomMatrixType & M, double value){
  /* add value to diagonal elements */
}

void update(ACustomStateType & v,        const double a,
		    const ACustomStateType & v1, const double b)
{
  // elementwise compute : v = a*v + b*v1
}

void update(ACustomStateType & v,        const double a,
            const ACustomStateType & v0, const double b,
            const ACustomStateType & v1, const double c)
{
  // elementwise compute : v = a*v + b*v0 + c*v1
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
#include "pressio/ode_steppers_implicit.hpp"

int main()
{}
@endcode

Obviously, if you want to use pressio nonlinear solvers, then you need provide
also the specializations described [here](md_pages_components_nonlinsolvers.html).

## API for discrete-time systems

@code{.cpp}
template<int num_states, class StateType, class SystemType>
auto create_arbitrary_stepper(const StateType & state,
						      SystemType && system);
@endcode

### Parameters

- `StateType`:
  - type of the data structure you use for the state

- `SystemType`:
  - class defining how to create and compute the residual and Jacobian

### Requirements

- `StateType`: must be copy constructible

- `SystemType` must conform to the following API:
  @code{.cpp}
  class ValidDiscreteTimeSystem
  {
	using scalar_type = /* whatever you need */;
	using state_type  = /* your type */;
	using discrete_time_residual_type = /* your type */;
	using discrete_time_jacobian_type = /* your type */;

	discrete_time_residual_type createDiscreteTimeResidual() const;
	discrete_time_jacobian_type createDiscreteTimeJacobian() const;

    // overload accepting 1 auxiliary state
	template<class StepCountType>
	void discreteTimeResidual(StepCountType,
							  scalar_type time,
							  scalar_type dt,
							  discrete_time_residual_type &,
							  const state_type &) const;

    // overload accepting 2 auxiliary states
	template<class StepCountType>
	void discreteTimeResidual(StepCountType,
								scalar_type time,
								scalar_type dt,
								discrete_time_residual_type &,
								const state_type &,
								const state_type &) const;

    // overload accepting 1 auxiliary state
	template<class StepCountType>
	void discreteTimeJacobian(StepCountType,
								scalar_type time,
								scalar_type dt,
								discrete_time_jacobian_type &,
								const state_type &) const;

    // overload accepting 2 auxiliary states
	template<class StepCountType>
	void discreteTimeJacobian(StepCountType,
								scalar_type time,
								scalar_type dt,
								discrete_time_jacobian_type &,
								const state_type &,
								const state_type &) const;
  };
  @endcode
