
#ifndef ROM_GALERKIN_UNSTEADY_DEFAULT_PROBLEM_MIXED_FOM_HPP_
#define ROM_GALERKIN_UNSTEADY_DEFAULT_PROBLEM_MIXED_FOM_HPP_

#include "mixed_fom_rom_unsteady_problem.hpp"

namespace pressio{ namespace rom{ namespace impl{

template <
  class IndVarType,
  class FomSystemType,
  class TrialSubspaceType,
  class GalerkinSystemType,
  template<class ...> class MixedFomRomStepper
  >
class GalerkinUnsteadyDefaultProblemRomFom
{

public:
  using independent_variable_type = IndVarType;
  using fom_state_type    = typename FomSystemType::state_type;
  using rom_state_type    = typename GalerkinSystemType::state_type;
  using rom_residual_type = typename GalerkinSystemType::rhs_type;
  using rom_jacobian_type = typename GalerkinSystemType::jacobian_type;

  using fom_stepper_type = decltype(ode::create_implicit_stepper(std::declval<ode::StepScheme>(),
   								 std::declval<const FomSystemType&>()));
  using rom_stepper_type = decltype(ode::create_implicit_stepper(std::declval<ode::StepScheme>(),
  								 std::declval<GalerkinSystemType>()));

  using mixed_stepper_t = MixedFomRomStepper<
    FomSystemType, fom_stepper_type, TrialSubspaceType, rom_stepper_type>;

  GalerkinUnsteadyDefaultProblemRomFom(ode::StepScheme odeSchemeName,
				       const TrialSubspaceType & trialSubspace,
				       const FomSystemType & fomSystem)
    : fomStepper_(ode::create_implicit_stepper(odeSchemeName, fomSystem))
    , trialSubspace_(trialSubspace)
    , romStepper_(ode::create_implicit_stepper(odeSchemeName, GalerkinSystemType(trialSubspace, fomSystem)))
    , mixed_(odeSchemeName, fomSystem, trialSubspace, fomStepper_, romStepper_)
  {
    assert(odeSchemeName_ == ode::StepScheme::BDF1 ||
  	   odeSchemeName_ == ode::StepScheme::BDF2);
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
  rom_stepper_type romStepper_;
  mixed_stepper_t mixed_;
};

}}} // end pressio::rom::impl
#endif  // ROM_IMPL_LSPG_UNSTEADY_PROBLEM_HPP_
