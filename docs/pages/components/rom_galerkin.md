
# rom: Galerkin: General Info

\todo: finish this, add more details

The pressio Galerkin ROMs are designed such that
using them involves these steps:

## 1. Create

You create an instance of one of the supported "Galerkin problem" using the following API:

```cpp
auto problem = pressio::rom::galerkin::create_<keywords>( /* args */ );
```

where `<keywords>` express the variant of the problem you want (more below),
and `args` are the arguments needed which depend on the variant you choose.
We currently support the following variants:

- Default: [link](md_pages_components_rom_galerkin_default.html)
- Hyper-reduced: [link](md_pages_components_rom_galerkin_hypred.html)
- Masked: [link](md_pages_components_rom_galerkin_masked.html)

The above `create` function returns a problem object that behaves like a stepper.
Therefore, you can use the problem like
you would with any other stepper object (more on this below).

### Explicit Problem

The problem meets the following API:

```cpp
class GalerkinExplicitProblemClass
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

where the traits contain:

```cpp
typename traits::fom_system_type;
typename traits::scalar_type;
typename traits::decoder_type;
typename traits::decoder_jac_type;
typename traits::galerkin_state_type;
typename traits::galerkin_velocity_type;
typename traits::stepper_type;
```

### Implicit Problem

The problem meets the following API:

```cpp
class GalerkinImplicitProblemClass
{
public:
  using traits = /* nested typedef to access the problem's traits */;

  // const ref to the object knowing how to reconstruct a FOM state
  const auto & fomStateReconstructor() const;

  using scalar_type    = /* ... */
  using state_type     = /* ... */;
  using residual_type  = /* ... */;
  using jacobian_type  = /* ... */;

  residual_type createResidual() const;
  jacobian_type createJacobian() const;
  void residual(const state_type& x, residual_type & res) const;
  void jacobian(const state_type& x, jacobian_type & jac) const;

  template <class TimeType, class StepCount, class ...Args>
  void operator()(galerkin_state_type & state,
				  const TimeType current_time,
				  const TimeType time_step_size_to_use,
				  const StepCount step_count,
				  Args && ... args);
};
```

where the traits contain:

```cpp
typename traits::fom_system_type;
typename traits::scalar_type;
typename traits::decoder_type;
typename traits::decoder_jac_type;
typename traits::galerkin_state_type;
typename traits::galerkin_residual_type;
typename traits::galerkin_jacobian_type;
typename traits::stepper_type;
```

## 2. Solve in time


### Example for explicit Galerkin

```cpp
int main()
{
  //...

  namespace pode = pressio::ode;
  namespace pgal = pressio::rom::galerkin;

  const auto scheme = pdoe::StepScheme:ForwardEuler;
  auto problem      = pgal::create_default_explicit_problem(scheme, /* args */);
  pode::advance_n_steps_and_observe(problem, /* args */);
}
```

### Example for implicit Galerkin

```cpp
int main()
{
  //...

  namespace pode = pressio::ode;
  namespace pgal = pressio::rom::galerkin;

  const auto scheme = pdoe::StepScheme:BDF1;
  auto problem      = pgal::create_default_implicit_problem(scheme, /* args */);
  auto solver = pressio::nonlinearsolvers::create_newton_rapshon(problem, /* args */);
  pode::advance_n_steps_and_observe(problem, /* args */, solver);
}
```


## Why does the problem behave like a stepper?

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


<!-- ## 2. Extract the stepper and solve in time -->

<!-- ```cpp -->
<!-- auto problem = pressio::rom::galerkin::create_<keyword>( /* args */ ); -->

<!-- // for explicit Galerkin, you can simply do for example: -->
<!-- pressio::ode::advance_n_steps_and_observe(stepper, ...); -->

<!-- // for implicit Galerking, you need to create a solver and solve -->
<!-- auto solver = pressio::solvers::create_newton_raphson(stepper, ...); -->

<!-- pressio::ode::advance_n_steps_and_observe(stepper, ..., solver); -->
<!-- ``` -->
