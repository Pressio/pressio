
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

The above `create` function returns a problem object that meets the following interface:

```cpp
class GalerkinProblemClass
{
public:
  using traits = /* nested typedef to access the problem's traits */;

  // returns the underlying stepper to use to solve the problem
  auto & stepper();

  // const ref to the object knowing how to reconstruct a FOM state
  const auto & fomStateReconstructor() const;
};
```

where the traits class contains the following:

```cpp
// for the explicit case, one can access the following traits:
typename traits::fom_system_type;
typename traits::scalar_type;
typename traits::decoder_type;
typename traits::decoder_jac_type;
typename traits::galerkin_state_type;
typename traits::galerkin_velocity_type;
typename traits::stepper_type;

// for the implicit case one has:
typename traits::fom_system_type;
typename traits::scalar_type;
typename traits::decoder_type;
typename traits::decoder_jac_type;
typename traits::galerkin_state_type;
typename traits::galerkin_residual_type;
typename traits::galerkin_jacobian_type;
typename traits::stepper_type;
```

The `stepper` method is, practically, what you would use
to retrieve the stepper and then use it to solve the problem.
The stepper method returns a non-const reference to an
[explicit stepper](md_pages_components_ode_steppers_explicit.html)
if, when you create the problem, you select an explicit scheme,
or an [implicit stepper](md_pages_components_ode_steppers_implicit.html)
if you select an implicit scheme.
Once you reference the stepper, you can then use it like
you would with any other stepper object (more on this below).

What does a stepper have to do with a Galerkin ROM?
The answer is that practically speaking, at the lowest-level,
a Galerkin problem can be reduced to simply a "custom" stepper to advance in time.
This is how pressio implements this and the reason why a Galerkin
problem contains a stepper object inside: when you create the
problem, pressio creates the appropriate custom stepper
object that you can use. You don't need to know how this is done,
or rely on the details, because these are problem- and implementation-dependent,
and we reserve the right to change this in the future.


## 2. Extract and Solve

### Example for explicit Galerkin

```cpp
int main()
{
  //...

  namespace pode = pressio::ode;
  namespace pgal = pressio::rom::galerkin;

  const auto scheme = pdoe::StepScheme:ForwardEuler;
  auto problem      = pgal::create_default_explicit_problem(scheme, /* args */);
  auto & stepper    = problem.stepper();

  pressio::ode::advance_n_steps_and_observe(stepper, /* args */);
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
  auto & stepper    = problem.stepper();

  auto solver = pressio::nonlinearsolvers::create_newton_rapshon(stepper, /* args */);
  pressio::ode::advance_n_steps_and_observe(stepper, /* args */, solver);
}
```


## Why not making the problem itself a stepper?

One could argue that the Galerkin problem class itself can be
made to behave like a stepper allowing one to avoid having to reference anything.
If we did this, then one would do, for example:

```cpp
auto problem = pgal::create_default_explicit_problem(scheme, /* args */);

// this is not valid, just for reasoning
pressio::ode::advance_n_steps_and_observe(problem, /* args */);
```

This is a valid comment, and one that we thought about.
We reserve the right to change this in the future, but for now
we have decided not to do this. We want to keep explicit
the need to reference the stepper once the problem is created.
This also makes it clear how to use the stepper.
This is just a choice, but we belive that keeping explicit the
dependence of Galerkin on time is important, and to clearly
convey that the problem "owns" something that you use.





<!-- ## 2. Extract the stepper and solve in time -->

<!-- ```cpp -->
<!-- auto problem = pressio::rom::galerkin::create_<keyword>( /* args */ ); -->

<!-- // for explicit Galerkin, you can simply do for example: -->
<!-- pressio::ode::advance_n_steps_and_observe(stepper, ...); -->

<!-- // for implicit Galerking, you need to create a solver and solve -->
<!-- auto solver = pressio::solvers::create_newton_raphson(stepper, ...); -->

<!-- pressio::ode::advance_n_steps_and_observe(stepper, ..., solver); -->
<!-- ``` -->
