
# rom: Unsteady LSPG

\todo: finish this, add more details

The pressio Unsteady LSPG ROMs are designed such that
using them involves these steps:

## 1. Create

You create an instance of one of the supported "LSPG problem" using the following API:

```cpp
auto problem = pressio::rom::lspg::create_<keywords>( /* args */ );
```

where `<keywords>` express the variant of the problem you want (more below),
and `args` are the arguments needed which depend on the variant you choose.
We currently support the following variants:

- Default: [link](md_pages_components_rom_lspg_default.html)
- Hyper-reduced: [link](md_pages_components_rom_lspg_hypred.html)
- Masked: [link](md_pages_components_rom_lspg_masked.html)

The above `create` function returns a problem object that meets the following interface:

```cpp
class UnsteadyLSPGProblemClass
{
public:
  using traits = /* nested typedef to access the problem's traits */;

  // returns the underlying stepper to use to solve the problem
  auto & stepper();

  // const ref to the object knowing how to reconstruct a FOM state
  const auto & fomStateReconstructor() const;
};
```

The `stepper` method is, practically, what you would use
to retrieve the stepper and then use it to solve the problem.
The stepper method returns a non-const reference to an
[implicit stepper](md_pages_components_ode_steppers_implicit.html).
Once you reference the stepper, you can then use it to solve the problem (more on this below).

What does a stepper have to do with a LSPG ROM?
The answer is that practically speaking, at the lowest-level,
a LSPG problem can be reduced to simply a "custom" stepper to advance in time.
When you create the
problem, pressio creates the appropriate custom stepper
object that you can use. You don't need to know how this is done,
or rely on the details, because these are problem- and implementation-dependent,
and we reserve the right to change this in the future.


## 2. Extract and Solve

```cpp
int main()
{
  //...

  namespace pode  = pressio::ode;
  namespace plspg = pressio::rom::lspg;

  const auto scheme = pdoe::StepScheme:BDF1;
  auto problem      = plspg::create_default_problem(scheme, /* args */);
  auto & stepper    = problem.stepper();

  auto solver = pressio::nonlinearsolvers::create_gauss_newton(stepper, /* args */);
  pressio::ode::advance_n_steps_and_observe(stepper, /* args */, solver);
}
```

Note that above we used as an example a Gauss-Newton solver because to solve
LSPG consists of "advancing" in time by solving at each step a
nonlinear least-squares problem.
