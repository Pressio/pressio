
// this is here so that we have a single place where
// we put the syntax and use it easily in other places
// without duplicating things

struct SyntaxOnly
{
  using independent_variable_type = /* your type */;
  using state_type                = /* your type */;
  using right_hand_side_type      = /* your type */;

  state_type createState() const;
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

  state_type createState() const;
  right_hand_side_type createRightHandSide() const;
  jacobian_type createJacobian() const;

  void rightHandSide(const state_type &,
                     const independent_variable_type &,
                     right_hand_side_type &) const;

  void jacobian(const state_type &,
                const independent_variable_type &,
                jacobian_type &) const;
};

struct SyntaxOnly
{
  using independent_variable_type = /* your type */;
  using state_type                = /* your type */;
  using right_hand_side_type      = /* your type */;
  using mass_matrix_type          = /* your type */;

  state_type createState() const;
  right_hand_side_type createRightHandSide() const;
  mass_matrix_type createMassMatrix() const;

  void rightHandSide(const state_type &,
                     const independent_variable_type &,
                     right_hand_side_type &) const;

  void massMatrix(const state_type &,
                  const independent_variable_type &,
                  mass_matrix_type &) const;
};

struct SyntaxOnly
{
  using independent_variable_type = /* your type */;
  using state_type                = /* your type */;
  using right_hand_side_type      = /* your type */;
  using mass_matrix_type          = /* your type */;

  state_type createState() const;
  right_hand_side_type createRightHandSide() const;
  mass_matrix_type createMassMatrix() const;

  void rightHandSide(const state_type &,
                     const independent_variable_type &,
                     right_hand_side_type &) const;

  void massMatrix(mass_matrix_type &) const;
};

struct SyntaxOnly
{
  using independent_variable_type = /* your type */;
  using state_type                = /* your type */;
  using right_hand_side_type      = /* your type */;
  using jacobian_type    = /* your type */;
  using mass_matrix_type = /* your type */;

  state_type createState() const;
  right_hand_side_type createRightHandSide() const;
  jacobian_type createJacobian() const;
  mass_matrix_type createMassMatrix() const;

  void rightHandSide(const state_type &,
                     const independent_variable_type &,
                     right_hand_side_type &) const;

  void jacobian(const state_type &,
                const independent_variable_type &,
                jacobian_type &) const;

  void massMatrix(const state_type &,
                  const independent_variable_type &,
                  mass_matrix_type &) const;
};

struct SyntaxOnly
{
  using independent_variable_type = /* your type */;
  using state_type                = /* your type */;
  using right_hand_side_type      = /* your type */;
  using jacobian_type    = /* your type */;
  using mass_matrix_type = /* your type */;

  state_type createState() const;
  right_hand_side_type createRightHandSide() const;
  jacobian_type createJacobian() const;
  mass_matrix_type createMassMatrix() const;

  void rightHandSide(const state_type &,
                     const independent_variable_type &,
                     right_hand_side_type &) const;

  void jacobian(const state_type &,
                const independent_variable_type &,
                jacobian_type &) const;

  void massMatrix(mass_matrix_type &) const;
};

struct FullyDiscreteSystem
{
  using independent_variable_type = /* ... */;
  using state_type                = /* ... */;
  using discrete_residual_type    = /* ... */;
  using discrete_jacobian_type    = /* ... */;

  state_type createState() const;
  discrete_residual_type createDiscreteResidual() const;
  discrete_jacobian_type createDiscreteJacobian() const;

  // accepting 2 states
  template<class StepCountType>
  void discreteResidual(StepCountType,
                        const independent_variable_type & /**/,
                        const independent_variable_type & /**/,
                        discrete_residual_type & /**/,
			discrete_jacobian_type &,
			bool computeJacobian,
                        const state_type & /**/,
			const state_type & /**/) const;

  // accepting 3 states
  template<class StepCountType>
  void discreteResidual(StepCountType,
                        const independent_variable_type & /**/,
                        const independent_variable_type & /**/,
                        discrete_residual_type &,
			discrete_jacobian_type &,
			bool computeJacobian,
                        const state_type &,
                        const state_type &,
                        const state_type &) const;

  // accepting 4 states
  template<class StepCountType>
  void discreteResidual(StepCountType,
                        const independent_variable_type & /**/,
                        const independent_variable_type & /**/,
                        discrete_residual_type &,
			discrete_jacobian_type &,
			bool computeJacobian,
                        const state_type &,
			const state_type &,
                        const state_type &,
                        const state_type &) const;
};






















class SyntaxOnly
{
  public:
    using state_type = /* some type */;
    using independent_variable_type  = /* some type */;

    void operator()(state_type & state,
                    pressio::ode::StepStartAt<independent_variable_type> /**/,
                    pressio::ode::StepCount /**/,
                    pressio::ode::StepSize<independent_variable_type> /**/);
};

class SyntaxOnly
{
  public:
    using state_type = /* some type */;
    using independent_variable_type  = /* some type */;

    template<class AuxT, class ...Args>
    void operator()(state_type & state,
                    pressio::ode::StepStartAt<independent_variable_type> /**/,
                    pressio::ode::StepCount /**/,
                    pressio::ode::StepSize<independent_variable_type> /**/,
                    AuxT && arg1,
                    Args && ... args);
};

class SyntaxOnly
{
  public:
    template<class IndepVarType>
    void operator()(const pressio::ode::StepCount /**/,
                    const pressio::ode::StepStartAt<IndepVarType> /**/,
                    pressio::ode::StepSize<IndepVarType> & /**/) const;
};

class SyntaxOnly
{
  public:
    template<class IndepVarType>
    void operator()(pressio::ode::StepCount /**/,
                    pressio::ode::StepStartAt<IndepVarType> /**/,
                    pressio::ode::StepSize<IndepVarType> & /**/,
                    pressio::ode::StepSizeMin<IndepVarType> & /**/t,
                    pressio::ode::StepSizeReduction<IndepVarType> & /**/) const;
};

class SyntaxOnly
{
  public:
    template<class IndepVarType, class StateType>
    void operator()(pressio::ode::StepCount /**/,
                    IndepVarType /**/,
                    const StateType & /**/) const;
};

class SyntaxOnly
{
  public:
    template<class IndepVarType, class StateType>
    void operator()(pressio::ode::StepCount /**/,
                    pressio::ode::StepStartAt<IndepVarType> /**/,
                    StateType & /**/) const;
};

struct ResidualAndJacobianPolicySyntaxOnly
{
  using independent_variable_type = /* your type */;
  using state_type                = /* your type */;
  using residual_type             = /* your type */;
  using jacobian_type             = /* your type */;

  state_type createState();
  residual_type createResidual();
  jacobian_type createJacobian();

  template <class StencilStatesContainerType, class StencilRhsContainerType>
  void operator()(pressio::ode::StepScheme                                   /**/,
		  const StateType &		                             /**/,
		  const StencilStatesContainerType &                         /**/,
		  StencilRhsContainerType &                                  /**/,
		  const pressio::ode::StepEndAt<independent_variable_type> & /**/,
		  pressio::ode::StepCount                                    /**/,
		  const pressio::ode::StepSize<independent_variable_type> &  /**/,
		  residual_type &		                             /**/,
		  jacobian_type &		                             /**/) const;
};

















// space above is left on purpose in case we need to
// expand something above wihout affecting all files that
// use literalinclude
struct SyntaxOnly
{
  using independent_variable_type = /* your type */;
  using state_type                = /* your type */;
  using right_hand_side_type      = /* your type */;
  using mass_matrix_type          = /* your type */;

  state_type createState() const;
  right_hand_side_type createRightHandSide() const;
  mass_matrix_type createMassMatrix() const;

  void rightHandSideAndMassMatrix(const state_type &,
				  const independent_variable_type &,
				  right_hand_side_type &,
				  mass_matrix_type &) const;
};
