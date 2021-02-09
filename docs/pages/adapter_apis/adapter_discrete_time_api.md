
# Discrete-time FOM adapter API

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
