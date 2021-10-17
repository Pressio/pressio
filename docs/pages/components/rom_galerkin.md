
# rom: Galerkin: high-level details

@m_class{m-note m-info}

@parblock
This page explains the API for using the pressio Galerkin ROMs.
After reading this, you should understand what a "pressio Galerkin problem" is,
the variants we currently support, and how to use the problem after instantiating it.

If anything is unclear, and/or you have suggestions on how
to improve this page, [open an issue on github](https://github.com/Pressio/pressio/issues).
@endparblock

<br/>

## Everything starts with creating a problem!

The main entry point to use the pressio Galerkin ROMs is the problem class.
You create an instance of one of the supported "Galerkin problems" as:

@m_class{m-block m-primary}

@par
```cpp
namespace pgal    = pressio::rom::galerkin;
const auto scheme = pressio::ode::StepScheme::some_value;
auto problem      = pgal::create_<keywords>(scheme, /* args */ );
```
@endparblock

where `<keywords>` express the variant you want (more below),
`some_value` is an enum value to select the time stepping scheme (e.g.,
explicit Euler, RK4, BDF1, etc), and `args` are the arguments needed which depend on the variant you choose.

We currently offer the following variants:

@m_div{m-button m-success}
<a href="md_pages_components_rom_galerkin_default.html">
@m_div{m-medium}&ensp;&emsp;Default Problem&emsp; &ensp; @m_enddiv
@m_div{m-small} click to learn more @m_enddiv
</a> @m_enddiv

@m_div{m-button m-primary}
<a href="md_pages_components_rom_galerkin_hypred.html">
@m_div{m-medium}Hyper-reduced Problem @m_enddiv
@m_div{m-small} click to learn more @m_enddiv
</a> @m_enddiv

@m_div{m-button m-warning}
<a href="md_pages_components_rom_galerkin_masked.html">
@m_div{m-medium}&ensp;&emsp; Masked Problem&ensp;&emsp; @m_enddiv
@m_div{m-small} click to learn more @m_enddiv
</a> @m_enddiv


## Problem class API

A Galerkin problem exposes different methods
depending on whether you choose an explicit of implicit scheme.
Below we discuss both scenarios.

### Explicit case

If you create an explicit Galerkin problem, the problem exposes the following API:

```cpp
class ExplicitGalerkinProblemClass
{
public:
  using traits = /* nested typedef to access the problem's traits */;

  // const ref to the object knowing how to reconstruct a FOM state
  const auto & fomStateReconstructor() const;

  template <class TimeType, class StepCount>
  void operator()(galerkin_state_type & state,
				  const TimeType current_time,
				  const TimeType time_step_size_to_use,
				  const StepCount step_count);
};
```

@m_class{m-block m-success}

@par Main thing to remember:
An explicit Galerkin problem satisfies the [steppable concept](md_pages_components_ode_advance.html)
(specifically, behaves like an [explicit stepper](md_pages_components_ode_steppers_explicit.html)).
@endparblock

If needed, the traits class contains:

```cpp
typename traits::fom_system_type;
typename traits::scalar_type;
typename traits::decoder_type;
typename traits::decoder_jac_type;
typename traits::galerkin_state_type;
typename traits::galerkin_velocity_type;
```

### Implicit case

If you create an implicit Galerkin problem, the problem exposes the following API:

```cpp
class ImplicitGalerkinProblemClass
{
public:
  using traits = /* nested typedef to access the problem's traits */;

  using scalar_type    = typename traits::scalar_type;
  using state_type     = typename traits::galerkin_state_type;
  using residual_type  = typename traits::galerkin_residual_type;
  using jacobian_type  = typename traits::galerkin_jacobian_type;

  // const ref to the object knowing how to reconstruct a FOM state
  const auto & fomStateReconstructor() const;

  template <class TimeType, class StepCount, class SolverType>
  void operator()(galerkin_state_type & state,
				  const TimeType current_time,
				  const TimeType time_step_size_to_use,
				  const StepCount step_count,
				  SolverType & solver)
  {
    // NOTE: here we reveal a bit of how operator() works.
	// This is on purpose so that some things below will be more clear.

    // The operator() has 3 main stages:
    // 1. we prepare for doing a step by setting/updating what is needed

    // 2. call solver (which will query operators)
	solver.solve(*this, state);

	// 3. we do some other things to end the step
  }

  residual_type createResidual() const;
  jacobian_type createJacobian() const;
  void residual(const state_type& x, residual_type & res) const;
  void jacobian(const state_type& x, jacobian_type & jac) const;
};
```

@m_class{m-block m-success}

@par Main thing to remember:
An implicit Galerkin problem satisfies the [steppable concept](md_pages_components_ode_advance.html)
(specifically, behaves like an [implicit stepper](md_pages_components_ode_steppers_implicit.html)).
@endparblock


The traits contain:

```cpp
typename traits::fom_system_type;
typename traits::scalar_type;
typename traits::decoder_type;
typename traits::decoder_jac_type;
typename traits::galerkin_state_type;
typename traits::galerkin_residual_type;
typename traits::galerkin_jacobian_type;
```

<br/>
___
<br/>

## How do I use an EXPLICIT problem?

See the following:

```cpp
namespace pode = pressio::ode;
namespace pgal = pressio::rom::galerkin;

const auto scheme = pdoe::StepScheme:ForwardEuler;
auto problem      = pgal::create_default_explicit_problem(scheme, /* args */);

// we assume a galState exists

// you can take a single step of the problem
const double currentTime = 0.;
const double dt = 0.5;
problem(galState, currentTime, dt, 1);

// or you can define your own stepping loop
double currentTime = 0.;
const double dt = 0.1;
for (int step = 1; step <= 10; step++){
 problem(galState, currentTime, dt, 1);
 currentTime += dt;
}

// or you can use our own functions to advance the problem
pode::advance_n_steps(problem, /* args */);

// or
pode::advance_n_steps_and_observe(problem, /* args */);

// or others
```

<br/>

## How do I use an IMPLICIT problem?

Recall that doing implicit time stepping it is not as simple as explicit.
[For implicit, in fact, you also need a *solver* to compute the solution at the next step](md_pages_components_ode_steppers_implicit.html).
In the case of Galerkin, you can use a Newton-Raphson solver,
because at eaach step, you are solving a (reduced) system
of equations with as many equations as the number of unknowns.
More specifically, the system you need to solve has as many equations as the
dimensionality of your approximating subspace.
See some sample snippets below:

```cpp
namespace psolvers = pressio::nonlinearsolvers;
namespace pode     = pressio::ode;
namespace pgal     = pressio::rom::galerkin;

const auto scheme = pdoe::StepScheme:BDF1;
auto problem      = pgal::create_default_implicit_problem(scheme, /* args */);

// we assume a galState exists

// for example, using the pressio solver
auto solver = psolvers::create_newton_rapshon(problem, /* args */);

// you can take a single step of the problem
{
  const double currentTime = 1.5;
  const double dt = 0.5;
  problem(galState, currentTime, dt, 1, solver);
}

// or you can define your own stepping loop
{
  double currentTime = 0.;
  const double dt = 0.1;
  for (int step = 1; step <= 10; step++){
   problem(galState, currentTime, dt, 1, solver);
   currentTime += dt;
  }
}

// or you can use our own functions to advance the problem
pode::advance_n_steps_and_observe(problem, /* args */, solver);
```

### Using a custom solver

In the snippet above, we show how to use the pressio
Newton-Raphson solver to solve the Galerkin problem.
If you want to use your own solver, you can do that!
Here we discuss how.

```cpp
template<class R_type, class J_type>
class CustomSolver
{
private:
  R_type m_R;
  J_type m_J;

public:
  template<class ProblemType>
  CustomSolver(ProblemType & system)
    : m_R(system.createResidual()),
	  m_J(system.createJacobian())
  {}

  template<class ProblemType, class StateType>
  void solve(ProblemType & problem, StateType & state)
  {
    // the problem object is the instance of the Galerkin problem you chose.
	// You can compute residual, Jacobian, etc and solve as you wish
	// for instance, you can do:

    for (/* some loop or whatever your solver needs */)
	{
	  problem.residual(state, m_R);
	  problem.jacobian(state, m_J);
  	  // do something with the residual (m_R) and Jacobian (m_J)
	}
  }
};

namespace pode     = pressio::ode;
namespace pgal     = pressio::rom::galerkin;

const auto scheme = pdoe::StepScheme:BDF1;
auto problem      = pgal::create_default_implicit_problem(scheme, /* args */);
using traits      = typename decltype(problem)::traits:

// instantiate your own solver
// and you can use it like shown above
using gal_R_type  = typename traits::galerkin_residual_type;
using gal_J_type  = typename traits::galerkin_jacobian_type;
CustomSolver<gal_R_type, gal_J_type> yourSolver(problem);

// we assume a galState exists

// you can take a single step of the problem
{
  const double currentTime = 0.;
  const double dt = 0.5;
  problem(galState, currentTime, dt, 1, yourSolver);
}

// or you can define your own stepping loop
{
  double currentTime = 0.;
  const double dt = 0.1;
  for (int step = 1; step <= 10; step++){
   problem(galState, currentTime, dt, 1, yourSolver);
   currentTime += dt;
  }
}

// or you can use our own functions to advance the problem
pode::advance_n_steps_and_observe(problem, /* args */, yourSolver);
```

@m_class{m-block m-warning}

@par Pay attention to:
If you need fine-grained access, i.e., you need to make the problem
take a single step (in other words you don't use our own `advance` methods),
you should *always* just rely on `operator()`.
To be more explicit, look at this snippet:
```cpp
// ...
// assume we have a galState
auto problem      = pgal::create_default_implicit_problem(scheme, /* args */);
solver = // create solver somehow

// !!! You should NOT do this, this is undefined behavior !!!
solver.solve(problem, galState)

// You should do this:
problem(..., solver);
// or, better, rely on (the advance functions do the right thing behind the scenes)
pode::advance_n_steps_and_observe(problem, /* args */, solver);
```
@endparblock


<br/>
___
<br/>

## Why does a Galerkin problem behave like a stepper?

The answer is that practically speaking, at the lowest-level,
a Galerkin problem can be reduced to simply a "custom" stepper to advance in time.
This is how pressio implements this and the reason why a Galerkin
problem contains a stepper object inside: when you create the
problem, pressio creates the appropriate custom stepper
object that you can use. You don't need to know how this is done,
or rely on the details, because these are problem- and implementation-dependent,
and we reserve the right to change this in the future.
All you need to know is that a an explicit Galerkin problem behaves like
a an explicit stepper, and an implicit Galerkin problem behaves
like an implicit stepper.
