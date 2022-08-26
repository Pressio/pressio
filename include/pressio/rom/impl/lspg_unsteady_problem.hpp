
#ifndef PRESSIO_ROM_IMPL_LSPG_PROBLEM_HPP_
#define PRESSIO_ROM_IMPL_LSPG_PROBLEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  class TrialSpaceType,
  class ResidualJacobianPolicyType
  >
class LspgUnsteadyProblem
{
  // note: to deduce the stepper_type it does not really matter
  // what scheme enum value we use, as long as it is an implicit one
  using stepper_type =
    decltype(::pressio::ode::create_implicit_stepper
	     (::pressio::ode::StepScheme::BDF1,
	      std::declval<ResidualJacobianPolicyType &>()
	      ));

  using fom_states_manager_type = LspgFomStatesManager<TrialSpaceType>;

public:
  using independent_variable_type  = typename ResidualJacobianPolicyType::independent_variable_type;
  using state_type    = typename TrialSpaceType::reduced_state_type;
  using residual_type = typename ResidualJacobianPolicyType::residual_type;
  using jacobian_type = typename ResidualJacobianPolicyType::jacobian_type;

  LspgUnsteadyProblem() = delete;

  template<class FomSystemType, class ...Args>
  LspgUnsteadyProblem(::pressio::ode::StepScheme odeSchemeName,
		      const TrialSpaceType & trialSpace,
		      const FomSystemType & fomSystem,
		      Args && ... args)
    : trialSpace_(trialSpace),
      fomStatesManager_(create_lspg_fom_states_manager(odeSchemeName, trialSpace)),
      rjPolicy_(trialSpace, fomSystem, fomStatesManager_, std::forward<Args>(args)...),
      stepper_( ::pressio::ode::create_implicit_stepper(odeSchemeName, rjPolicy_))
  {}

  template<class SolverType, class ...ArgsOp>
  void operator()(state_type & reducedState,
		  pressio::ode::StepStartAt<independent_variable_type> sStart,
		  pressio::ode::StepCount sCount,
		  pressio::ode::StepSize<independent_variable_type> sSize,
		  SolverType & solver,
		  ArgsOp && ...argsop)
  {
    stepper_(reducedState, sStart, sCount, sSize,
     	     solver, std::forward<ArgsOp>(argsop)...);
  }

  auto createState() const{
    return stepper_.createState();
  }

  residual_type createResidual() const{
    return stepper_.createResidual();
  }

  jacobian_type createJacobian() const{
    return stepper_.createJacobian();
  }

  void residualAndJacobian(const state_type & odeState,
			   residual_type & R,
			   jacobian_type & J,
			   bool computeJacobian) const{
    stepper_.residualAndJacobian(odeState, R, J, computeJacobian);
  }

private:
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  LspgFomStatesManager<TrialSpaceType> fomStatesManager_;
  ResidualJacobianPolicyType rjPolicy_;
  stepper_type stepper_;
};

}}} // end pressio::rom::impl
#endif
