
#ifndef PRESSIO_ROM_IMPL_LSPG_PROBLEM_HPP_
#define PRESSIO_ROM_IMPL_LSPG_PROBLEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  class TrialSpaceType,
  class ResidualPolicyType,
  class JacobianPolicyType
  >
class LspgUnsteadyProblem
{
  // note: to deduce the stepper_type it does not really matter
  // what scheme enum value we use, as long as it is an implicit one
  using stepper_type =
    decltype(::pressio::ode::create_implicit_stepper
	     (::pressio::ode::StepScheme::BDF1,
	      std::declval<ResidualPolicyType &>(),
	      std::declval<JacobianPolicyType &>()
	      ));

  using fom_states_manager_type = LspgFomStatesManager<TrialSpaceType>;

public:
  using independent_variable_type  = typename ResidualPolicyType::independent_variable_type;
  using state_type    = typename TrialSpaceType::reduced_state_type;
  using residual_type = typename ResidualPolicyType::residual_type;
  using jacobian_type = typename JacobianPolicyType::jacobian_type;

  LspgUnsteadyProblem() = delete;

  template<class FomSystemType, class ...Args>
  LspgUnsteadyProblem(::pressio::ode::StepScheme odeSchemeName,
		      const TrialSpaceType & trialSpace,
		      const FomSystemType & fomSystem,
		      Args && ... args)
    : trialSpace_(trialSpace),
      fomStatesManager_(create_lspg_fom_states_manager(odeSchemeName, trialSpace)),
      resPolicy_(trialSpace, fomSystem, fomStatesManager_, std::forward<Args>(args)...),
      jacPolicy_(trialSpace, fomSystem, fomStatesManager_, std::forward<Args>(args)...),
      stepper_( ::pressio::ode::create_implicit_stepper(odeSchemeName, resPolicy_, jacPolicy_))
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

  void residual(const state_type & odeState, residual_type & R) const{
    stepper_.residual(odeState, R);
  }

  void jacobian(const state_type & odeState, jacobian_type & J) const{
    stepper_.jacobian(odeState, J);
  }

private:
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  LspgFomStatesManager<TrialSpaceType> fomStatesManager_;
  ResidualPolicyType resPolicy_;
  JacobianPolicyType jacPolicy_;
  stepper_type stepper_;
};

}}} // end pressio::rom::impl
#endif
