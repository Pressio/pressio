
// this is here so that we have a single place where
// we put the syntax and use it easily in other places
// without duplicating things

struct SyntaxOnly
{
  using time_type            = /* your type */;
  using state_type           = /* your type */;
  using right_hand_side_type = /* your type */;

  right_hand_side_type createRightHandSide() const;

  void rightHandSide(const state_type &,
		     const time_type &,
		     right_hand_side_type &) const;
};


struct SyntaxOnly
{
  using time_type            = /* your type */;
  using state_type           = /* your type */;
  using right_hand_side_type = /* your type */;

  right_hand_side_type createRightHandSide() const;

  void rightHandSide(const state_type &,
                     const time_type &,
                     right_hand_side_type &) const;

  template<class BasisType>
  /*return type*/ createApplyJacobianResult(const BasisType & ) const;

  template<class BasisType, class ResultType>
  void applyJacobian(const state_type &,
		     const BasisType &,
		     const time_type &,
		     ResultType &) const;
};


class SyntaxOnly
{
  public:
    using state_type    = /* your type */;
    using residual_type = /* your type */;

  public:
    residual_type createResidual() const;

    template<class BasisType>
    /*return type*/ createApplyJacobianResult(const BasisType & B) const;

    void residual(const state_type &,
		  residual_type &) const;

    template<class BasisType, class ResultType>
    void applyJacobian(const state_type & state,
		       const BasisType &,
		       ResultType &) const;
};


struct SyntaxOnly
{
  using time_type              = /* ... */;
  using state_type             = /* ... */;
  using discrete_residual_type = /* ... */;

  discrete_residual_type createDiscreteTimeResidual() const;

  template<class BasisType>
  /* return type */ createResultOfDiscreteTimeJacobianActionOn(const BasisType &) const;

  // accepting 2 states
  template<class StepCountType>
  void discreteTimeResidualAndJacobianAction(StepCountType /*stepCount*/,
					     const time_type & /*t_n+1*/,
					     const time_type & /*dt*/,
					     discrete_residual_type & /*R*/,
					     const /* basis type */ &,
					     bool computeJacobianAction,
					     /* jacobian action result's type */ &,
					     const state_type & /*y_n+1*/,
					     const state_type & /*y_n*/) const;

  // accepting 3 states
  template<class StepCountType>
  void discreteTimeResidualAndJacobianAction(StepCountType /*stepCount*/,
					     const time_type & /*t_n+1*/,
					     const time_type & /*dt*/,
					     discrete_residual_type & /*R*/,
					     const /* basis type */ &,
					     bool computeJacobianAction,
					     /* jacobian action result's type */ &,
					     const state_type & /*y_n+1*/,
					     const state_type & /*y_n*/,
					     const state_type & /*y_n-1*/) const;
};
