
# Adapter/system class API for Unsteady LSPG

For unsteady, Pressio supports a so-called *continuous-time* and *discrete-time* API. \todo 

## Continuous-time API

The continuous-time API operates such that the user is responsible 
to compute the continuous-time operators, e.g., the velocity, and pressio assembles the 
discrete-time operators. It is an API that very expressive of the formulation 
on which pressio relies on.

### Basic Unsteady LSPG 

For a basic unsteady LSPG with the *continuous-time* API we have:

@m_class{m-code-figure} @parblock
@code{.cpp}
class AdapterUnsteadyLSPG{
public:
  using scalar_type =
  using state_type =
  using velocity_type =
  using dense_matrix_type =

public:
  velocity_type createVelocity() const;
  dense_matrix_type createApplyJacobianResult(const dense_matrix_type &) const;

  void velocity(state, time, velo) const;
  void applyJacobian(state, const dense_matrix_type & B, time, dense_matrix_type & A) const; // computes: A = Jac B
};
@endcode

### Preconditioned LSPG
If you want a preconditioned problem, then the above class must be extended to add: 
@m_class{m-code-figure} @parblock
@code{.cpp}
class AdapterSteadyLSPG
{
  // everything from above

  // for preconditioned problem
  void applyPreconditioner(const state_type&, const scalar_type &t, velocity_type & maskedObj) const;
  void applyPreconditioner(const state_type&, const scalar_type &t, dense_matrix_type & discreteTimeJ) const;

};
@endcode
<!--   // for preconditioned problem
  // for masked problem
  residual_type createApplyMaskResult(const residual_type & unmaskedObj) const;
  dense_matrix_type createApplyMaskResult(const dense_matrix_type & unmaskedObj) const;
  void applyMask(const unmaskedObj, residual_type & maskedObj) const;
  void applyMask(const unmaskedObj, dense_matrix_type & maskedObj) const;
 -->


## Discrete-time API

@m_class{m-code-figure} @parblock
@code{.cpp}
class 
{
public:
  using scalar_type = //..;
  using state_type  = //...;
  using discrete_time_residual_type = //...;
  using dense_matrix_type = //...;

public:
  discrete_time_residual_type createDiscreteTimeResidual() const;
  dense_matrix_type createApplyDiscreteTimeJacobianResult(const dense_matrix_type &) const
  { // let A =  tdJac * B
    dense_matrix_type A(/* construct A */);
    return A;
  }

  template <typename step_t, typename ... Args>
  void discreteTimeResidual(const step_t & step,
            const scalar_type & time,
            const scalar_type & dt,
            discrete_time_residual_type & R,
            pressio::Norm normKind,
            scalar_type & normR,
	          //variadic # of states (user sets stencil size)
            Args & ... states) const
  {
    this->discreteTimeResidualImpl(step, time, dt, R, std::forward<Args>(states)... );
  }

  template <typename step_t, typename ... Args>
  void applyDiscreteTimeJacobian(const step_t & step,
           const scalar_type & time,
           const scalar_type & dt,
           const dense_matrix_type & B,
           dense_matrix_type & A,
	         //variadic # of states (user sets stencil size)
           Args & ... states) const
  {
    this->applyDiscreteTimeJacobianImpl(step, time, dt, B, stateIdForJacobian,
					A, std::forward<Args>(states)...);
  }
};
@endcode

