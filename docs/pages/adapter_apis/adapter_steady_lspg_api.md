
# FOM Adapter API for Steady LSPG

Note that the adapter classes shown below serve the purpose of interfacing
your native application with pressio, but the actual object instantiated
from one of these should be created as part of your application.
These classes do **not** contain anything strictly related to pressio,
but just contain types native to your application.
\todo (fix)

## Basic Steady LSPG problem

The adapter class should look like:

@m_class{m-code-figure} @parblock
@code{.cpp}
class AdapterSteadyLSPG
{
  // ...
public:
  // The following aliases MUST be exposed because Pressio detects them.
  // If these are not visible, mispelled or not found, you get a compile-time error
  // because your adapter class does not the right API
  using scalar_type       = /* your native scalar type */
  using state_type        = /* your native state type */
  using residual_type     = /* your native residual type */
  using dense_matrix_type = /* your native dense matrix type */

public:
  // creates the residual object
  // This is only called once to create the operators, does not need to contain real data.
  residual_type createResidual() const;

  // creates the result of applying the jacobian to the argument
  // This is only called once to create the operators, does not need to contain real data.
  dense_matrix_type createApplyJacobianResult(const dense_matrix_type &) const;

  void residual(state, r) const;

  // computes the result of applying the jacobian to the argument: A  = Jacobian B
  void applyJacobian(state, B, A) const; // computes: A = Jac B
};
@endcode
@endparblock

<br>
### Preconditioned steady LSPG

If you want a preconditioned steady LSPG problem, then the above class must be extended to add:
@m_class{m-code-figure} @parblock
@code{.cpp}
class AdapterSteadyLSPG
{
  // everything from above

  void applyPreconditioner(const state_type &, residual_type & r) const;
  void applyPreconditioner(const state_type &, dense_matrix_type & jac) const;
};
@endcode
@endparblock

<!--   // for preconditioned problem
  // for masked problem
  residual_type createApplyMaskResult(const residual_type & unmaskedObj) const;
  dense_matrix_type createApplyMaskResult(const dense_matrix_type & unmaskedObj) const;
  void applyMask(const unmaskedObj, residual_type & maskedObj) const;
  void applyMask(const unmaskedObj, dense_matrix_type & maskedObj) const;
 -->
