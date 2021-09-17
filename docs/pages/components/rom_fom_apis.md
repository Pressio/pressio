
# FOM adapter API

In order to use the `presssio/rom` model reduction functionalities,
obviously there needs to be a way to exchange
data/information between pressio and your FOM application.
To do so, pressio requires you to write an *adapter class* as
a minimally intrusive layer to standardize the way pressio interfaces with any application.
As explained below, in general, preparing the adapter should only
involve *exposing* from your applications some operators.

## Stead API

Intended for when your FOM application is expressed as
@f[
\boldsymbol{R}(\boldsymbol{y}; \boldsymbol{\mu}) = 0
@f]
where @f$y@f$ is the full-order model (FOM) state,
@f$R@f$ is the residual
\todo finish.

### Synopsis

```cpp
class FomSteadyAdapter
{
public:
  using scalar_type   = /* your type */;
  using state_type    = /* your type */;
  using residual_type = /* your type */;

public:
  // you need to create an instance of the residual
  residual_type createResidual() const;

  // given a state, you compute the residual
  void residual(const state_type &,
				residual_type &) const;
};
```

### Usage
- @m_span{m-text m-warning}This adapter can ONLY be used for doing stedy LSPG ROMs.@m_endspan
- See the following examples: \toadd

<br/>


## Continuous-time API: Velocity Only

Intended for when your FOM application is expressed in *time-continuous* form as
@f[
\frac{d \boldsymbol{y}}{dt} =
\boldsymbol{f}(\boldsymbol{y},t; \boldsymbol{\mu}),
\quad \boldsymbol{y}(0;\boldsymbol{\mu}) = \boldsymbol{y}(\boldsymbol{\mu}),
@f]
where @f$y@f$ is the full-order model (FOM) state,
@f$f@f$ is what we call the FOM velocity (or RHS), and @f$t@f$ is time,
and, for some reason, you can/want to only expose
the right-hand-side (or velocity) of your FOM application.
\todo finish.

### Synopsis

```cpp
class FomAdapter
{
public:
  using scalar_type   = /* your type */;
  using state_type    = /* your type */;
  using velocity_type = /* your type */;

public:
  // you need to create an instance of the velocity (or RHS)
  velocity_type createVelocity() const;

  // given a state and time, you compute the velocity (or RHS)
  void velocity(const state_type &,
			    const scalar_type & time,
				velocity_type &) const;
};
```

### Usage
- @m_span{m-text m-warning}This adapter can ONLY be used for doing Galerkin ROMs with explicit time stepping.@m_endspan
- See the following examples: ...

<br/>

## Continuous-time API: Velocity and Jacobian action

This API is intended for any system expressible in *time-continuous* form as above,
but you expose both the right-hand-side of your FOM application as well as
the action of the velocity's Jacobian on some operand (more on this later).

### Synopsis
```cpp
class FomAdapter
{
public:
  using scalar_type   = /* your type */;
  using state_type    = /* your type */;
  using velocity_type = /* your type */;

public:
  velocity_type createVelocity() const;

  void velocity(const state_type &,
			    const scalar_type & time,
				velocity_type &) const;

  // OperandType: typically the data type you use for the decoder's jacobian
  template<class OperandType>
  OperandType createApplyJacobianResult(const OperandType &) const;

  // computes: A = Jac B
  template<class OperandType>
  void applyJacobian(const state_type &,
					 const OperandType & B,
					 const scalar_type & time,
					 OperandType & A) const;
};
```

### Usage

- @m_span{m-text m-warning}can be used for doing Galerkin ROMs with explicit and implicit time stepping@m_endspan,
- @m_span{m-text m-warning}can be used for LSPG and WLS (note that LSPG and WLS only
make sense for implicit time integration).@m_endspan
- See the following examples: ...

<br/>



## Discrete-time API

This API is intended for any system expressible in a discrete-time form as
@f[
\boldsymbol{R}(\boldsymbol{y_{n+1}}, \boldsymbol{y_{n}}, \boldsymbol{y_{n-1}}, ..., t_{n+1}, dt_{n+1}; ...) = \boldsymbol{0}
@f]
where @f$y@f$ is the full-order model (FOM) state, @f$t@f$ is time, and @f$R@f$ is the residual.
\todo finish.

### Synopsis

```cpp
class
{
public:
  using scalar_type					= /* your type */;
  using state_type				    = /* your type */;
  using discrete_time_residual_type = /* your type */;

public:
  discrete_time_residual_type createDiscreteTimeResidual() const;

  // OperandType: typically the data type you use for the decoder's jacobian
  template<class OperandType>
  operand_t createApplyDiscreteTimeJacobianResult(const OperandType &) const;

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

  // OperandType: typically the data type you use for the decoder's jacobian
  template <class step_t, OperandType, class ... Args>
  void applyDiscreteTimeJacobian(const step_t & step,
								 const scalar_type & time,
								 const scalar_type & dt,
								 const OperandType & B,
								 OperandType & A,
								 const state_type & y_np1,
								 const state_type & y_n) const
  {
    // given y_n+1, y_n
	// compute A = dR/dy_n+1 B
  }

  // OperandType: typically the data type you use for the decoder's jacobian
  template <typename step_t, typename ... Args>
  void applyDiscreteTimeJacobian(const step_t & step,
								 const scalar_type & time,
								 const scalar_type & dt,
								 const OperandType & B,
								 OperandType & A,
								 const state_type & y_np1,
								 const state_type & y_n,
								 const state_type & y_nm1) const
  {
    // given y_n+1, y_n, y_n-1
	// compute A = dR/dy_n+1 B
  }
};
```

### Usage

- @m_span{m-text m-warning}for doing Galerkin *implicit* time stepping@m_endspan
- @m_span{m-text m-warning}for doing LSPG @m_endspan
- See the following examples: ...

<br/>


## Should you prefer the continuous-time or discrete-time API?

In general, we suggest users to always prefer the continuous-time API because it is more general.
However, there are situations where the discrete-time API is more useful or even necessary.
For example, \todo.
