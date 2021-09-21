
# Nonlinear Solvers - Problem APIs


## Residual-Jacobian API

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

## Hessian-Gradient API
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
