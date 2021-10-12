
# Nonlinear Solvers: Levenberg-Marquardt


@m_class{m-note m-default}

@parblock
Defined in header `<pressio/solvers_nonlinear.hpp>`

Public namespace: `pressio::nonlinearsolvers`
@endparblock


## Levenberg–Marquardt

### API, Parameters and Requirements

```cpp
template<
  class ProblemClassType,
  class StateType,
  class LinearSolverType
  >																(1)
auto create_levenberg_marquardt(const ProblemClassType & system,
                                const StateType & state,
                                LinearSolverType && lsolver);

template<
  class ProblemClassType,
  class StateType,
  class LinearSolverType,
  class WeightingOpType
  >																(2)
auto create_levenberg_marquardt(const ProblemClassType & system,
                                const StateType & state,
                                LinearSolverType && lsolver,
						        WeightingOpType && weightOperator);
```

- `system`:
  - instance of the problem you want to solve
	@m_class{m-block m-warning}

	@parblock
	- overload 1: accepts the [residual-jacobian, hessian-gradient API, or their fused versions](md_pages_components_nonlinsolvers_system_api.html)
	- overload 2: *only* accepts the [residual-jacobian API or its fused version](md_pages_components_nonlinsolvers_system_api.html)
	@endparblock

- `state`:
  - your state data
  - Requirements: must be an Eigen or Kokkos vector: \todo explain why

- `lsolver`:
  - linear solver for solving the normal equations, choose one from [linear solver API](md_pages_components_linsolvers.html)
  - if you want to implement your own, then the linear solver class still has to conform to the [linear solver API](md_pages_components_linsolvers.html)

- `weightOperator`:
  - weighting operator for doing weighted least-squares.
  Must conform to:
  ```cpp
  class WeightingOperator
  {
    public:
 	void operator()(const residual_type & operand, residual_type & result);
	void operator()(const jacobian_type & operand, jacobian_type & result);
  };
  ```


### Ops
### Example usage
