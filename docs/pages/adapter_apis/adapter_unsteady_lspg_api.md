
# FOM Adapter API for Unsteady LSPG

\todo say that for LSPG, the adapater is the one for galerkin plus something

## API for Basic Unsteady LSPG Problem

For a basic unsteady LSPG with the *continuous-time* API we have:

@m_class{m-code-figure} @parblock
@code{.cpp}
class AdapterUnsteadyLSPG{
public:
  using scalar_type =
  using state_type =
  using velocity_type =

public:
  velocity_type createVelocity() const;
  void velocity(state, time, velo) const;

  // operand_type should be the data (matrix) type you used to store the basis.
  operand_t createApplyJacobianResult(const operand_t &) const;

  // computes: A = Jac B
  void applyJacobian(state, const operand_t & B, time, operand_t & A) const;
};
@endcode
@endparblock

<!-- ## Preconditioned LSPG -->
<!-- If you want a preconditioned problem, then the above class must be extended to add: -->
<!-- @m_class{m-code-figure} @parblock -->
<!-- @code{.cpp} -->
<!-- class AdapterSteadyLSPG -->
<!-- { -->
<!--   // everything from above -->

<!--   // for preconditioned problem -->
<!--   void applyPreconditioner(const state_type&, const scalar_type &t, velocity_type & maskedObj) const; -->
<!--   void applyPreconditioner(const state_type&, const scalar_type &t, dense_matrix_type & discreteTimeJ) const; -->

<!-- }; -->
<!-- @endcode -->
<!-- @endparblock -->

<!-- <\!--   // for preconditioned problem -->
<!--   // for masked problem -->
<!--   residual_type createApplyMaskResult(const residual_type & unmaskedObj) const; -->
<!--   dense_matrix_type createApplyMaskResult(const dense_matrix_type & unmaskedObj) const; -->
<!--   void applyMask(const unmaskedObj, residual_type & maskedObj) const; -->
<!--   void applyMask(const unmaskedObj, dense_matrix_type & maskedObj) const; -->
<!--  -\-> -->

## Discrete-time API

The discrete-time API for LSPG is the same as for all methods,
and you can see it [here](./md_pages_adapter_apis_adapter_discrete_time_api.html).
