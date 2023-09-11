
#ifndef ROM_MIXED_FOM_ROM_UNSTEADY_PROBLEM_HPP_
#define ROM_MIXED_FOM_ROM_UNSTEADY_PROBLEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  class FomSystemType,
  class FomStepperType,
  class TrialSubspaceType,
  class RomStepperType
  >
class MixedFomRomStepper
{
  enum class StepKind{ fom, rom, none };

public:
  static_assert( std::is_same_v<
		 typename FomStepperType::independent_variable_type,
		 typename RomStepperType::independent_variable_type>);
  using independent_variable_type = typename FomStepperType::independent_variable_type;
  using fom_state_type            = typename FomStepperType::state_type;
  using rom_state_type            = typename RomStepperType::state_type;

  MixedFomRomStepper(ode::StepScheme odeSchemeName,
		     const FomSystemType & fomSystem,
		     const TrialSubspaceType & trialSubspace,
		     FomStepperType & fomStepper,
		     RomStepperType & romStepper)
    : odeSchemeName_(odeSchemeName)
    , fomStepper_(fomStepper)
    , fomAuxStateForTransitionFomToRom_(fomSystem.createState())
    , fomAuxStateForTransitionRomToFom_(fomSystem.createState())
    , fomScratchState_(fomSystem.createState())
    , trialSubspace_(trialSubspace)
    , romStepper_(romStepper)
    , romAuxStateForTransitionRomToFom_(trialSubspace.createReducedState())
  {
    assert(odeSchemeName_ == ode::StepScheme::BDF1 ||
  	   odeSchemeName_ == ode::StepScheme::BDF2);
  }

  template<class SolverType, class ...ArgsOp>
  void operator()(impl::FomStepTag /*unused*/,
		  fom_state_type & fomState,
		  pressio::ode::StepStartAt<independent_variable_type> sStart,
		  pressio::ode::StepCount sCount,
		  pressio::ode::StepSize<independent_variable_type> sSize,
		  SolverType & solver,
		  ArgsOp && ...argsop)
  {

    //preconditions
    if (activeFlag_ == StepKind::rom){
      const std::string msg = "You cannot take a FOM step starting from a ROM one. " \
	"You must first transition from ROM to FOM, and then you can execute a regular FOM step";
      throw std::runtime_error(msg);
    }

    PRESSIOLOG_DEBUG("FOM is active: doing **FOM** step = ", sCount.get(), "\n");
    activeFlag_ = StepKind::fom;

    if (odeSchemeName_ == ode::StepScheme::BDF2){
      pressio::ops::deep_copy(fomAuxStateForTransitionFomToRom_, fomState);
    }
    fomStepper_(fomState, sStart, sCount, sSize,
		solver, std::forward<ArgsOp>(argsop)...);
  }

  template<class SolverType, class ...ArgsOp>
  void operator()(impl::RomStepTag /*unused*/,
		  rom_state_type & romState,
		  pressio::ode::StepStartAt<independent_variable_type> sStart,
		  pressio::ode::StepCount sCount,
		  pressio::ode::StepSize<independent_variable_type> sSize,
		  SolverType & solver,
		  ArgsOp && ...argsop)
  {

    //preconditions
    if (activeFlag_ == StepKind::fom){
      const std::string msg = "You cannot take a ROM step starting from a FOM one. "
	"You must first transition from FOM to ROM,  and then you can execute a regular ROM step";
      throw std::runtime_error(msg);
    }

    PRESSIOLOG_DEBUG("ROM is active: doing **ROM** step = ", sCount.get(), "\n");
    activeFlag_ = StepKind::rom;

    if (odeSchemeName_ == ode::StepScheme::BDF2){
      pressio::ops::deep_copy(romAuxStateForTransitionRomToFom_, romState);
    }
    romStepper_(romState, sStart, sCount, sSize,
		solver, std::forward<ArgsOp>(argsop)...);
  }

  template<class SolverType, class ...ArgsOp>
  void operator()(impl::TransitionToRomAndDoStepTag /*unused*/,
		  fom_state_type & fomState,
		  rom_state_type & romState,
		  pressio::ode::StepStartAt<independent_variable_type> sStart,
		  pressio::ode::StepCount sCount,
		  pressio::ode::StepSize<independent_variable_type> sSize,
		  SolverType & solver,
		  ArgsOp && ...argsop)
  {

    //preconditions
    if (activeFlag_ != StepKind::fom){
      throw std::runtime_error("Transitioning from FOM to ROM must be done from a FOM context.");
    }

    PRESSIOLOG_DEBUG("Trasitioning FOM->ROM, and doing **ROM** step = ", sCount.get(), "\n");
    activeFlag_ = StepKind::rom;

    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    const auto & shift = trialSubspace_.get().translationVector();

    // 1. project the FOM state to compute the rom state to start from
    // romState = phi^T (fomState - shift)
    pressio::ops::update(fomScratchState_, 0, fomState, 1, shift, -1);
    pressio::ops::product(::pressio::transpose(), 1., phi, fomScratchState_, 0., romState);

    // 2. take the ROM step
    if (odeSchemeName_ == ode::StepScheme::BDF1){
      // things are easy, I don't need to overwrite any history since BDF1 only depends on current state
      romStepper_(romState, sStart, sCount, sSize, solver, std::forward<ArgsOp>(argsop)...);
    }
    else{
      assert(odeSchemeName_ == ode::StepScheme::BDF2);

      // I have to save the state for transitioning that might be used later
      pressio::ops::deep_copy(romAuxStateForTransitionRomToFom_, romState);

      // I have to modify the ROM state stored in the ode stencils at n-1 by doing
      auto FomToRomTransition = [&](auto & odeStencilStates){
	auto romStateTmp = pressio::ops::clone(romState);
	pressio::ops::update(fomScratchState_, 0, fomAuxStateForTransitionFomToRom_, 1, shift, -1);
	pressio::ops::product(::pressio::transpose(), 1., phi, fomScratchState_, 0., romStateTmp);

	auto & y_nm1 = odeStencilStates(ode::nMinusOne());
	pressio::ops::deep_copy(y_nm1, romStateTmp);
      };

      // call the stepper
      romStepper_(romState, sStart, sCount, sSize, FomToRomTransition,
		  solver, std::forward<ArgsOp>(argsop)...);
    }
  }

  template<class SolverType, class ...ArgsOp>
  void operator()(impl::TransitionToFomAndDoStepTag /*unused*/,
		  fom_state_type & fomState,
		  rom_state_type & romState,
		  pressio::ode::StepStartAt<independent_variable_type> sStart,
		  pressio::ode::StepCount sCount,
		  pressio::ode::StepSize<independent_variable_type> sSize,
		  SolverType & solver,
		  ArgsOp && ...argsop)
  {

    //preconditions
    if (activeFlag_ != StepKind::rom){
      throw std::runtime_error("Transitioning from ROM to FOM must be done from a ROM context.");
    }

    PRESSIOLOG_DEBUG("Trasitioning ROM->FOM, and doing **FOM** step = ", sCount.get(), "\n");
    activeFlag_ = StepKind::fom;

    // 1. reconstruct the fom state from the current rom state
    trialSubspace_.get().mapFromReducedState(romState, fomState);

    // 2. take the FOM step
    if (odeSchemeName_ == ode::StepScheme::BDF1){
      // things are easy, I don't need to overwrite any history since BDF1 only depends on current state
      fomStepper_(fomState, sStart, sCount, sSize, solver, std::forward<ArgsOp>(argsop)...);
    }
    else{
      assert(odeSchemeName_ == ode::StepScheme::BDF2);

      // I have to save the state for transitioning that might be used later
      pressio::ops::deep_copy(fomAuxStateForTransitionFomToRom_, fomState);

      // I have to modify the state stored at n-1 using the projection of the FOM
      auto RomToFomTransition = [&](auto & odeStencilStates){
	trialSubspace_.get().mapFromReducedState( romAuxStateForTransitionRomToFom_, fomScratchState_);
	auto & y_nm1 = odeStencilStates(ode::nMinusOne());
	pressio::ops::deep_copy(y_nm1, fomScratchState_);
      };

      // call the stepper
      fomStepper_(fomState, sStart, sCount, sSize, RomToFomTransition,
		  solver, std::forward<ArgsOp>(argsop)...);
    }
  }

private:
  StepKind activeFlag_ = StepKind::none;
  ode::StepScheme odeSchemeName_;

  // members for fom
  std::reference_wrapper<FomStepperType> fomStepper_;
  fom_state_type fomAuxStateForTransitionFomToRom_;
  fom_state_type fomAuxStateForTransitionRomToFom_;
  fom_state_type fomScratchState_;

  // members for rom
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<RomStepperType> romStepper_;
  rom_state_type romAuxStateForTransitionRomToFom_;
};

}}} // end pressio::rom::impl
#endif  // ROM_IMPL_LSPG_UNSTEADY_PROBLEM_HPP_
