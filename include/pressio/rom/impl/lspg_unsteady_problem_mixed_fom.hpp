
#ifndef ROM_IMPL_LSPG_UNSTEADY_PROBLEM_FOMROM_HPP_
#define ROM_IMPL_LSPG_UNSTEADY_PROBLEM_FOMROM_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  class IndepVarType,
  class FomSystemType,
  class TrialSubspaceType,
  class ResJacPolicyType,
  template<class ...> class MixedFomRomStepper
  >
class LspgUnsteadyProblemRomFom
{

public:
  using independent_variable_type = IndepVarType;
  using fom_state_type   = typename FomSystemType::state_type;
  using rom_state_type   = typename TrialSubspaceType::reduced_state_type;
  using fom_stepper_type = decltype(ode::create_implicit_stepper(ode::StepScheme::BDF1,
								 std::declval<const FomSystemType&>()));
  using rom_stepper_type = decltype(ode::create_implicit_stepper(ode::StepScheme::BDF1,
								 std::declval<ResJacPolicyType &>()));

  using mixed_stepper_t = MixedFomRomStepper<
    FomSystemType, fom_stepper_type, TrialSubspaceType, rom_stepper_type>;

  template<class ...Args>
  LspgUnsteadyProblemRomFom(::pressio::ode::StepScheme odeSchemeName,
			    const TrialSubspaceType & trialSubspace,
			    const FomSystemType & fomSystem,
			    Args && ... args)
    : fomStepper_(ode::create_implicit_stepper(odeSchemeName, fomSystem))
    , trialSubspace_(trialSubspace)
    , fomStatesManager_(create_lspg_fom_states_manager(odeSchemeName, trialSubspace))
    , romStepperPolicy_(trialSubspace, fomSystem, fomStatesManager_, std::forward<Args>(args)...)
    , romStepper_(ode::create_implicit_stepper(odeSchemeName, romStepperPolicy_))
    , mixed_(odeSchemeName, fomSystem, trialSubspace, fomStepper_, romStepper_)
  {
    assert(odeSchemeName == ode::StepScheme::BDF1 ||
	   odeSchemeName == ode::StepScheme::BDF2);
  }

  template<class Tag, class StateType, class SolverType, class ...ArgsOp>
  mpl::enable_if_t< std::is_same_v<Tag, FomStepTag> || std::is_same_v<Tag, RomStepTag> >
  operator()(Tag tag,
	     StateType & fomOrRomState,
	     pressio::ode::StepStartAt<independent_variable_type> sStart,
	     pressio::ode::StepCount sCount,
	     pressio::ode::StepSize<independent_variable_type> sSize,
	     SolverType & solver,
	     ArgsOp && ...argsop)
  {
    static_assert(std::is_same_v<StateType, rom_state_type> ||
		  std::is_same_v<StateType, fom_state_type>);
    mixed_(tag, fomOrRomState, sStart, sCount, sSize,
	   solver, std::forward<ArgsOp>(argsop)...);
  }

  template<class Tag, class SolverType, class ...ArgsOp>
  mpl::enable_if_t< std::is_same_v<Tag, TransitionToRomAndDoStepTag> ||
		    std::is_same_v<Tag, TransitionToFomAndDoStepTag> >
  operator()(Tag tag,
	     fom_state_type & fomState,
	     rom_state_type & romState,
	     pressio::ode::StepStartAt<independent_variable_type> sStart,
	     pressio::ode::StepCount sCount,
	     pressio::ode::StepSize<independent_variable_type> sSize,
	     SolverType & solver,
	     ArgsOp && ...argsop)
  {
    mixed_(tag, fomState, romState, sStart, sCount, sSize,
    	   solver, std::forward<ArgsOp>(argsop)...);
  }

private:
  fom_stepper_type fomStepper_;
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  LspgFomStatesManager<TrialSubspaceType> fomStatesManager_;
  ResJacPolicyType romStepperPolicy_;
  rom_stepper_type romStepper_;
  mixed_stepper_t mixed_;
};

}}} // end pressio::rom::impl
#endif  // ROM_IMPL_LSPG_UNSTEADY_PROBLEM_HPP_
