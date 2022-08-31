
// this is here so that we have a single place where
// we put the syntax and use it easily in other places
// without duplicating things

struct SyntaxOnly
{
  using independent_variable_type = /* your type */;
  using state_type                = /* your type */;
  using right_hand_side_type      = /* your type */;

  state_type           createState() const;
  right_hand_side_type createRightHandSide() const;

  void rightHandSide(const state_type &,
		     const independent_variable_type &,
		     right_hand_side_type &) const;
};

struct SyntaxOnly
{
  using independent_variable_type = /* your type */;
  using state_type                = /* your type */;
  using right_hand_side_type      = /* your type */;
  using jacobian_type             = /* your type */;

  state_type	       createState() const;
  right_hand_side_type createRightHandSide() const;
  jacobian_type	       createJacobian() const;

  void rightHandSide(const state_type &,
                     const independent_variable_type &,
                     right_hand_side_type &) const;

  void jacobian(const state_type &,
                const independent_variable_type &,
                jacobian_type &) const;
};

struct FullyDiscreteSystem
{
  using step_count_type = typename pressio::ode::StepCount::value_type;
  using independent_variable_type = /* ... */;
  using state_type                = /* ... */;
  using discrete_residual_type    = /* ... */;
  using discrete_jacobian_type    = /* ... */;

  state_type             createState() const;
  discrete_residual_type createDiscreteResidual() const;
  discrete_jacobian_type createDiscreteJacobian() const;

  // 2 total states
  void discreteResidualAndJacobian(step_count_type /*stepNumber*/,
				   const independent_variable_type & /*t_n+1*/,
				   const independent_variable_type & /*dt*/,
				   discrete_residual_type & /*R*/,
				   discrete_jacobian_type & /*J*/,
				   bool computeJacobian,
				   const state_type & /*statePredictionAt_n+1*/,
				   const state_type & /*stateAt_n */) const;


  // 3 total states
  void discreteResidualAndJacobian(step_count_type /*stepNumber*/,
				   const independent_variable_type & /*t_n+1*/,
				   const independent_variable_type & /*dt*/,
				   discrete_residual_type & /*R*/,
				   discrete_jacobian_type & /*J*/,
				   bool computeJacobian,
				   const state_type & /*statePredictionAt_n+1*/,
				   const state_type & /*stateAt_n*/,
				   const state_type & /*stateAt_n-1*/) const;

  // 3 total states
  void discreteResidualAndJacobian(step_count_type /*stepNumber*/,
				   const independent_variable_type & /*t_n+1*/,
				   const independent_variable_type & /*dt*/,
				   discrete_residual_type & /*R*/,
				   discrete_jacobian_type & /*J*/,
				   bool computeJacobian,
				   const state_type & /*statePredictionAt_n+1*/,
				   const state_type & /*stateAt_n*/,
				   const state_type & /*stateAt_n-1*/,
				   const state_type & /*stateAt_n-2*/) const;
};
