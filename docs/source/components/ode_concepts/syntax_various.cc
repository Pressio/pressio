
// this is here so that we have a single place where
// we put the syntax and use it easily in other places
// without duplicating things

// steppable
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

// steppable with args
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

// step size policy
class SyntaxOnly
{
  public:
    template<class IndepVarType>
    void operator()(const pressio::ode::StepCount /**/,
                    const pressio::ode::StepStartAt<IndepVarType> /**/,
                    pressio::ode::StepSize<IndepVarType> & /**/) const;
};

// step size policy with reduction
class SyntaxOnly
{
  public:
    template<class IndepVarType>
    void operator()(pressio::ode::StepCount /**/,
                    pressio::ode::StepStartAt<IndepVarType> /**/,
                    pressio::ode::StepSize<IndepVarType> & /**/,
                    pressio::ode::StepSizeMin<IndepVarType> & /**/,
                    pressio::ode::StepSizeScalingFactor<IndepVarType> & /**/) const;
};

// state obeserver
class SyntaxOnly
{
  public:
    template<class IndepVarType, class StateType>
    void operator()(pressio::ode::StepCount /**/,
                    IndepVarType /**/,
                    const StateType & /**/) const;
};

// state guesser
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
