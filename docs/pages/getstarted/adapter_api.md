
# FOM Adapter API

In order to use the ROM functionalities in pressio, 
one needs to create a FOM adapter class.
Pressio supports two main types of API as described below.

## Continuous-time API

Recall that pressio targets a generic full-order model (FOM) system written as
@f[
\frac{d \boldsymbol{y}}{dt} =
\boldsymbol{f}(\boldsymbol{y},t; \boldsymbol{\mu}),
\quad \boldsymbol{y}(0;\boldsymbol{\mu}) = \boldsymbol{y}(\boldsymbol{\mu}),
@f]
where @f$y@f$ is the FOM state and @f$f(...)@f$ is the FOM velocity.

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
              //variadic # of states (user sets stencil size)
              Args & ... states) const
  {
    this->discreteTimeResidualImpl(step, time, dt, R, std::forward<Args>(states)... );
  }

  template <typename step_t, typename ... Args>
  void applyDiscreteTimeJacobian(const step_t & step,
                 const scalar_type & time,
                 const scalar_type & dt,
                 const operand_t & B,
                 operand_t & A,
                 //variadic # of states (user sets stencil size)
                 Args & ... states) const
  {
    this->applyDiscreteTimeJacobianImpl(step, time, dt, B, stateIdForJacobian,
          A, std::forward<Args>(states)...);
  }
};
```
