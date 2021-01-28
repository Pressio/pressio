
# Continuous-time FOM adapter API


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

## A: velocity only
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

@par Where can you do with `AdapterA`?
This kind of adapter can be used for doing Galerkin ROMs with explicit time stepping.


## B: velocity and Jacobian action
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

@par Where can you do with `AdapterB`?
This kind of adapter can be used for doing Galerkin ROMs with explicit
and implicit time stepping, LSPG and WLS (note that LSPG and WLS only
make sense for implicit time integration).
