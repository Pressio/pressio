
# FOM Adapter discrete-time API

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
@endparblock
