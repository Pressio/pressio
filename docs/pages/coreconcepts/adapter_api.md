
# Adapter API

@m_class{m-block m-info}

@par
An adapter class allows a FOM application to expose data via an API conforming to Pressio requirements.

<!-- To use the functionalities in pressio, obviously there needs to be
a way to exchange data/information between pressio and your FOM application.
To do so, in pressio we leverage the idea of an *adapter class* as a layer
allowing to standardize the way pressio interfaces with any application.
Schematically, the flow of interfation is shown below:
@image html schem.svg width=65%
 -->

## Continuous-time API

This API is intended for any system expressible in *time-continuous* form as
@f[
\frac{d \boldsymbol{y}}{dt} =
\boldsymbol{f}(\boldsymbol{y},t; \boldsymbol{\mu}),
\quad \boldsymbol{y}(0;\boldsymbol{\mu}) = \boldsymbol{y}(\boldsymbol{\mu}),
@f]
where @f$y@f$ is the full-order model (FOM) state,
@f$f@f$ the FOM velocity, and @f$t@f$ is time.
\todo finish.

We envision two scenarios:
* A: you are only able to expose the right-hand-side (or velocity) of your FOM application
* B: you expose the right-hand-side of your FOM application as well as
the action of the velocity's Jacobian on some operand

### A: velocity only
```cpp
class AdapterA
{
public:
  using scalar_type =
  using state_type =
  using velocity_type =

public:
  velocity_type createVelocity() const;

  void velocity(const state_type &,
			    const scalar_type & time,
				velocity_type &) const;
};
```

@m_class{m-block m-warning}

@par Where can you use the `AdapterA`?
This version of the adapter can be used for doing Galerkin ROMs with explicit time stepping.


### B: velocity and Jacobian action
```cpp
class AdapterB
{
public:
  using scalar_type =
  using state_type =
  using velocity_type =

public:
  velocity_type createVelocity() const;

  void velocity(const state_type &,
			    const scalar_type & time,
				velocity_type &) const;

  // operand_type is the data (matrix) type you use
  // to represent the decoder's jacobian
  operand_t createApplyJacobianResult(const operand_t &) const;

  // computes: A = Jac B
  void applyJacobian(const state_type &,
					 const operand_t & B,
					 const scalar_type & time,
					 operand_t & A) const;
};
```

@m_class{m-block m-warning}

@par Where can you use the `AdapterB`?
This version of the adapter can be used for doing Galerkin ROMs with explicit
and implicit time stepping, LSPG and WLS (note that LSPG and WLS only
make sense for implicit time integration).


## Discrete-time API

This API is intended for any system expressible in a discrete-time form as
@f[
\boldsymbol{R}(\boldsymbol{y}, \boldsymbol{y_{n-1}}, ..., t_n, dt_n; ...) = \boldsymbol{0}
@f]
where @f$y@f$ is the full-order model (FOM) state, @f$t@f$ is time, and @f$R@f$ is the residual.
\todo finish.

```cpp
class
{
public:
  using scalar_type = //..;
  using state_type  = //...;
  using discrete_time_residual_type = //...;

public:
  discrete_time_residual_type createDiscreteTimeResidual() const;

  // operand_type should be the data (matrix) type you used to store the basis.
  operand_t createApplyDiscreteTimeJacobianResult(const operand_t &) const
  {
    // let A =  discreteTimeJac * B
    operand_t A(/* construct A */);
    return A;
  }

  template <typename step_t, typename ... Args>
  void discreteTimeResidual(const step_t & step,
                            const scalar_type & time,
							const scalar_type & dt,
							discrete_time_residual_type & R,
							const state_type & y_np1,
							const state_type & y_n) const
  {
    // given y_n+1, y_n
	// compute R
  }

  template <typename step_t, typename ... Args>
  void discreteTimeResidual(const step_t & step,
                            const scalar_type & time,
							const scalar_type & dt,
							discrete_time_residual_type & R,
							const state_type & y_np1,
							const state_type & y_n,
							const state_type & y_nm1) const
  {
    // given y_n+1, y_n, y_n-1
	// compute R
  }


  template <typename step_t, typename ... Args>
  void applyDiscreteTimeJacobian(const step_t & step,
								 const scalar_type & time,
								 const scalar_type & dt,
								 const operand_t & B,
								 operand_t & A,
								 const state_type & y_np1,
								 const state_type & y_n) const
  {
    // given y_n+1, y_n
	// compute A = dR/dy_n+1 B
  }

  template <typename step_t, typename ... Args>
  void applyDiscreteTimeJacobian(const step_t & step,
								 const scalar_type & time,
								 const scalar_type & dt,
								 const operand_t & B,
								 operand_t & A,
								 const state_type & y_np1,
								 const state_type & y_n,
								 const state_type & y_nm1) const
  {
    // given y_n+1, y_n, y_n-1
	// compute A = dR/dy_n+1 B
  }
};
```

@m_class{m-block m-warning}

@par Where can you use the discrete-time API?
This version of the adapter can be **only** used for doing Galerkin and LSPG ROMs with *implicit* time stepping.



@m_class{m-block m-info}

@par Should one prefer the continuous-time or discrete-time API?
In general, we suggest users to always prefer the continuous-time API because it is more general.



<!--

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

public:
  // creates the residual object
  // This is only called once to create the operators, does not need to contain real data.
  residual_type createResidual() const;

  // creates the result of applying the jacobian to the operand.
  // This is only called once to create the operators, does not need to contain real data.
  // operand_type should be the data (matrix) type you used to store the basis.
  operand_type createApplyJacobianResult(const operand_type &) const;

  void residual(state, r) const;

  // computes the result of applying the jacobian to the argument: A  = Jacobian B
  void applyJacobian(state, B, A) const; // computes: A = Jac B
};
@endcode
@endparblock -->
