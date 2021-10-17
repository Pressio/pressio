
# rom: Steady LSPG

\todo: write this better

The pressio steady LSPG ROMs are designed around two main steps:

## 1. Create

You instantiate a "steady LSPG problem", e.g.:

@m_class{m-block m-primary}

@par
```cpp
namespace plspg = pressio::rom::lspg;
auto problem    = plspg::create_keyword_steady_problem(/* args */);
```
@endparblock

where `keyword` expresses the variant of the problem you want (more below).

We currently support three variants:

@m_div{m-button m-success}
<a href="md_pages_components_rom_lspg_default_steady.html">
@m_div{m-medium}&ensp;&emsp;Default Problem&emsp; &ensp; @m_enddiv
@m_div{m-small} click to learn more @m_enddiv
</a> @m_enddiv

@m_div{m-button m-primary}
<a href="md_pages_components_rom_lspg_hypred_steady.html">
@m_div{m-medium}Hyper-reduced Problem @m_enddiv
@m_div{m-small} click to learn more @m_enddiv
</a> @m_enddiv

@m_div{m-button m-warning}
<a href="md_pages_components_rom_lspg_masked_steady.html">
@m_div{m-medium}&ensp;&emsp; Masked Problem&ensp;&emsp; @m_enddiv
@m_div{m-small} click to learn more @m_enddiv
</a> @m_enddiv


Refer to each problem page for details on each specific variant.

The returned `problem` object is an instantiation of a class exposing the following interface:

```cpp
class SteadyLspgProblem
{
public:
  using traits = /* nested typedef with trait class */;

  using scalar_type    = /* ... */;
  using state_type     = /* ... */;
  using residual_type  = /* ... */;
  using jacobian_type  = /* ... */;

  residual_type createResidual() const;
  jacobian_type createJacobian() const;
  void residual(const state_type& x, residual_type & res) const;
  void jacobian(const state_type& x, jacobian_type & jac) const;

  // const ref to the object knowing how to reconstruct a FOM state
  const auto & fomStateReconstructor() const;
};
```

## 2. Solve

- you instantiate and use a nonlinear least-squares solver of your choice to solve the problem.
  Note, in fact, that the problem's API conforms to the one required by the nonlinear solvers

- for this solve stage, you don't have to use the pressio4py solvers.
  Once you have the problem object, you can also use your own nonlinear least-squares solver.
  As shown above, the `problem` exposes all the operators that you need to solve.
