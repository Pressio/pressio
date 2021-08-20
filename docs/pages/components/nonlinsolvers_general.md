
# Nonlinear Solvers - General

Defined in header `<pressio/solvers_nonlinear.hpp>`

Public namespace: `pressio::nonlinearsolvers`

## Overview

| Name                                  | Use for:                                          | Useful math links     |
|-------------------------------------- |---------------------------------------------  |---------------- |
| Newton-Raphson                        | Systems of nonlinear equations                | [l1](https://link.springer.com/content/pdf/bbm%3A978-3-319-69407-8%2F1.pdf), [l2](https://www.cmu.edu/math/undergrad/suami/pdfs/2014_newton_method.pdf)   |
| Gauss-Newton                          | Nonlinear least-squares problem.              | [l1](https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm)                |
| Levenberg–Marquardt                   | Nonlinear least-squares problem.              | [l1](https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm)                 |
| Iteratively reweighted least squares  | optimization problem formulated in a p-norm   | [l1](https://en.wikipedia.org/wiki/Iteratively_reweighted_least_squares) |


## General solver settings

- convergence criterion
- logging
- updating step, choices and custom


## Problem APIs

```cpp
struct ProblemResJac
{
  using scalar_type    = /* your type */
  using state_type     = /* your type */;
  using residual_type  = /* your type */;
  using jacobian_type  = /* your type */;

  residual_type createResidual() const;
  jacobian_type createJacobian() const;
  void residual(const state_type& x, residual_type & res) const;
  void jacobian(const state_type& x, jacobian_type & jac) const;
};
```

```cpp
struct ProblemFusedResJac
{
  using scalar_type    = /* your type */
  using state_type     = /* your type */;
  using residual_type  = /* your type */;
  using jacobian_type  = /* your type */;

  residual_type createResidual() const;
  jacobian_type createJacobian() const;
  void residualAndJacobian(const state_type& x,
					       residual_type & res,
						   jacobian_type & jac,
						   bool recomputeJacobian) const;
};
```

```cpp
struct ProblemHessGrad
{
  using scalar_type   = /* your type */
  using state_type    = /* your type */;
  using hessian_type  = /* your type */;
  using gradient_type = /* your type */;

  hessian_type createHessian() const;
  gradient_type createGradient() const;

  void residualNorm(const state_type & state,
					pressio::Norm normKind,
				    scalar_type & resNorm) const;

  void gradient(const state_type &,
				gradient_type &,
				pressio::Norm normKind,
				scalar_type & normResidual,
				bool recomputeJacobian) const;

  void hessian(const state_type &, hessian_type &) const;
};
```

```cpp
struct ProblemHessGradFused
{
  using scalar_type   = /* your type */
  using state_type    = /* your type */;
  using hessian_type  = /* your type */;
  using gradient_type = /* your type */;

  hessian_type createHessian() const;
  gradient_type createGradient() const;

  void residualNorm(const state_type & state,
					pressio::Norm normKind,
				    scalar_type & resNorm) const;

  void hessianAndGradient(const state_type &,
						  hessian_type &,
						  gradient_type &,
						  pressio::Norm normKind,
						  scalar_type & normResidual,
						  bool recomputeJacobian) const;
};
```
