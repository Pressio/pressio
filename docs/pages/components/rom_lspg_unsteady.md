
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

The above `create` function returns a problem object that behaves like a stepper.
Therefore, you can use the problem like
you would with any other stepper object (more on this below).


```cpp
class UnsteadyLSPGProblemClass
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
  void operator()(lspg_state_type & state,
				  const TimeType current_time,
				  const TimeType time_step_size_to_use,
				  const StepCount step_count,
				  Args && ... args);
};
```

## 2. Solve in time

```cpp
int main()
{
  //...

  namespace pode  = pressio::ode;
  namespace plspg = pressio::rom::lspg;

  const auto scheme = pdoe::StepScheme:BDF1;
  auto problem      = plspg::create_default_problem(scheme, /* args */);
  auto solver = pressio::nonlinearsolvers::create_gauss_newton(stepper, /* args */);
  pode::advance_n_steps_and_observe(problem, /* args */, solver);
}
```

Note that above we used as an example a Gauss-Newton solver because to solve
LSPG consists of "advancing" in time by solving at each step a
nonlinear least-squares problem.
