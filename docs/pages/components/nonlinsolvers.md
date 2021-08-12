
# Nonlinear Solvers

Defined in header `<pressio/solvers_nonlinear.hpp>`

Public namespace: `pressio::nonlinearsolvers`


## Overview

| Name                                	| Use for:                                         	| Useful math links   	|
|--------------------------------------	|---------------------------------------------	|----------------	|
| Newthon-Raphson                      	| Systems of nonlinear equations              	| [l1](https://link.springer.com/content/pdf/bbm%3A978-3-319-69407-8%2F1.pdf), [l2](https://www.cmu.edu/math/undergrad/suami/pdfs/2014_newton_method.pdf) 	|
| Gauss-Newton (normal equations)      	| Nonlinear least-squares problem.            	| [l1](https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm)               	|
| Gauss-Newton (via QR)                	| Nonlinear least-squares problem.            	| [l1](https://en.wikipedia.org/wiki/Gauss%E2%80%93Newton_algorithm)               	|
| Levenberg–Marquardt                  	| Nonlinear least-squares problem.            	| [l1](https://en.wikipedia.org/wiki/Levenberg%E2%80%93Marquardt_algorithm)               	|
| Iteratively reweighted least squares 	| optimization problem formulated in a p-norm 	| [l1](https://en.wikipedia.org/wiki/Iteratively_reweighted_least_squares) |

## Newthon-Raphson

### API, Parameters and Requirements
```cpp
template<class ProblemClassType, class StateType, class LinearSolverType>
auto create_newton_raphson(const ProblemClassType & system,
						   const StateType & state,
						   LinearSolverType && lsolver);
```

- `ProblemClassType`:
  - type of the class defining the problem to solve
  - Requirements: must conform to either the residual-jacobian API or the fused residual-jacobian API

- `StateType`:
  - type of the data structure you use for the state
  - Requirements: must be copy constructible

- `LinearSolverType`:
  - self-explanatory
  - Requirements: must conform to this [API](/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio/docs/html/md_pages_components_linsolvers.html)

### Ops
When using custom data types not supported in [pressio ops](/Users/fnrizzi/Desktop/work/ROM/gitrepos/pressio/docs/html/md_pages_components_ops.html), you need to specialize a trait class and some operations
such that pressio can operate on your data. For the sake of explanation, suppose that you use:

```cpp
using scalar_type   = double;
using state_type    = ACustomVectorClass;  //this can be any type
using jacobian_type = ACustomMatrixClass;  //this can be any type
```

Then you need to provide the following specializations:

```cpp
namespace pressio{

template<> struct Traits<ACustomVectorClass>{
  using scalar_type = double;
};

template<> struct Traits<ACustomMatrixClass>{
  using scalar_type = double;
};

namespace ops{

ACustomVectorClass clone(const ACustomVectorClass & src){ /* return a deep copy of src */ }
ACustomMatrixClass clone(const ACustomMatrixClass & src){ /* return a deep copy of src */ }

void set_zero(ACustomVectorClass & object){ /* set elements zero */ }
void set_zero(ACustomMatrixClass & object){ /* set elements zero */ }

scalar_type norm2(const ACustomVectorClass & object){
  /* return l2-norm of object */
}

void update(ACustomVectorClass & v,		   scalar_type a,
			const ACustomVectorClass & v1, scalar_type b)
{
  // compute v = a*v + v1*b;
}

void scale(ACustomVectorClass & v, scalar_type factor){
  /* scale v elementwise by factor
}

}}//end namepsace pressio::ops
```

### Example usage
```cpp
int main()
{
  // assuming that:
  // problem_t is a problem class that meets API
  // state_t is defined too

  problem_t myProblem;
  state_t y(10);
  // set initial state

  // create linear system
  using lin_solver_t = /* something that meets API needed */;
  lin_solver_t linearSolverObj;

  namespace pnls = pressio::nonlinearsolvers;
  auto NonLinSolver = pnls::create_newton_raphson(myProblem, y, linearSolverObj);
  NonLinSolver.solve(myProblem, y);
}
```

## Gauss-Newton via Normal-Equations

### API, Parameters and Requirements
### Ops
### Example usage


## Gauss-Newton via QR factorization

### API, Parameters and Requirements
### Ops
### Example usage



## Levenberg–Marquardt

### API, Parameters and Requirements
### Ops
### Example usage



## General solver settings

- convergence criterion
- logging
- updating step, choices and custom



## Problem Classes

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
