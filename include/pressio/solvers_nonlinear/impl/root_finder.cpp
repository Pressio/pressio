
#ifndef SOLVERS_NONLINEAR_IMPL_ROOT_FINDER_HPP_
#define SOLVERS_NONLINEAR_IMPL_ROOT_FINDER_HPP_

namespace pressio{
namespace nonlinearsolvers{
namespace impl{

template<
  class ProblemTag,
  class UserDefinedSystemType,
  class RegistryType,
  class ToleranceType,
  class NormDiagnosticsContainerType,
  class DiagnosticsLoggerType,
  class UpdaterType>
void root_solving_loop_impl(ProblemTag /*problemTag*/,
          const UserDefinedSystemType & system,
          RegistryType & reg,
          Stop stopEnumValue,
          ToleranceType stopTolerance,
          NormDiagnosticsContainerType & normDiagnostics,
          const DiagnosticsLoggerType & logger,
          int maxIters,
          UpdaterType && updater)
{

  using state_type = typename UserDefinedSystemType::state_type;

  auto objective = [&reg, &system](const state_type & stateIn){
    auto & r = reg.template get<ResidualTag>();
    system.residualAndJacobian(stateIn, r, {});
    return ::pressio::ops::norm2(r);
  };

  auto mustStop = [
		   &normDiag = std::as_const(normDiagnostics),
		   stopEnumValue, maxIters, stopTolerance]
    (int stepCount)
  {
    const Diagnostic stopDiag = stop_criterion_to_public_diagnostic(stopEnumValue);
    switch (stopEnumValue){
    case Stop::AfterMaxIters:
      return stepCount == maxIters;
    default:
      if (is_absolute_diagnostic(stopDiag)){
	return normDiag[stopDiag].getAbsolute() < stopTolerance;
      }else{
	return normDiag[stopDiag].getRelative() < stopTolerance;
      }
    };
  };

  int iStep = 0;
  while (++iStep <= maxIters){
    /* stage 1 */
    try{
      compute_residual_and_jacobian(reg, system);
    }
    catch (::pressio::eh::ResidualEvaluationFailureUnrecoverable const &e){
      PRESSIOLOG_CRITICAL(e.what());
      throw ::pressio::eh::NonlinearSolveFailure();
    }
    catch (::pressio::eh::ResidualHasNans const &e){
      PRESSIOLOG_CRITICAL(e.what());
      throw ::pressio::eh::NonlinearSolveFailure();
    }

    /* stage 2 */
    solve_newton_step(reg);

    /* stage 3 */
    std::for_each(normDiagnostics.begin(), normDiagnostics.end(),
      [&reg, iStep](auto & v){
        const bool isFirstIteration = iStep==1;
        compute_norm_internal_diagnostics(reg, isFirstIteration, v);
      });
    logger(iStep, normDiagnostics);

    /* stage 4*/
    if (mustStop(iStep)){
      PRESSIOLOG_DEBUG("nonlinsolver: stopping");
      break;
    }

    /* stage 5 */
    try{
      const auto currentObjValue =
  normDiagnostics[InternalDiagnostic::residualAbsoluteRelativel2Norm].getAbsolute();
      updater(reg, objective, currentObjValue);
    }
    catch (::pressio::eh::LineSearchStepTooSmall const &e) {
      // nicely exist the solve
      PRESSIOLOG_WARN(e.what());
      break;
    }
    catch (::pressio::eh::LineSearchObjFunctionChangeTooSmall const &e) {
      // nicely exist the solve
      PRESSIOLOG_WARN(e.what());
      break;
    }
  }
}


template<class Tag, class StateType, class RegistryType, class NormValueType>
class RootFinder : public RegistryType
{
  Tag tag_;
  int maxIters_ = 100;
  Stop stopEnValue_ = Stop::WhenAbsolutel2NormOfCorrectionBelowTolerance;
  NormValueType stopTolerance_ = 0.000001;
  Update updateEnValue_ = Update::Standard;

  using norm_diagnostics_container = DiagnosticsContainer<
    InternalDiagnosticDataWithAbsoluteRelativeTracking<NormValueType> >;
  norm_diagnostics_container normDiagnostics_;
  DiagnosticsLogger diagnosticsLogger_ = {};

public:
  template<class ...Args>
  RootFinder(Tag tagIn,
	     const std::vector<Diagnostic> & diags,
	     Args && ...args)
    : RegistryType(std::forward<Args>(args)...),
      tag_(tagIn),
      normDiagnostics_(diags)
  {

    // currently we don't have the diagonostics stuff all flushed out
    // so we limit it to work for a specific case
    const auto & publicDiags = normDiagnostics_.publicNames();
    assert(publicDiags.size() == 4);
    assert(publicDiags[0] == Diagnostic::residualAbsolutel2Norm);
    assert(publicDiags[1] == Diagnostic::residualRelativel2Norm);
    assert(publicDiags[2] == Diagnostic::correctionAbsolutel2Norm);
    assert(publicDiags[3] == Diagnostic::correctionRelativel2Norm);

    // check the stop criterion uses a metric already supported in the diagonstics
    const auto stopMetric = stop_criterion_to_public_diagnostic(stopEnValue_);
    normDiagnostics_.addIfUnsupported(stopMetric);
    // need to reset the logger since the names might have changed
    diagnosticsLogger_.resetFor(publicDiags);
  }

  // query/set update criterion
  Update currentUpdateCriterion() const   { return updateEnValue_; }
  void setUpdateCriterion(Update value) { updateEnValue_ = value; }

  // query/set stop criterion, tolerance
  Stop currentStopCriterion() const          { return stopEnValue_; }
  void setStopCriterion(Stop value)	     { stopEnValue_ = value; }
  void setStopTolerance(NormValueType value) { stopTolerance_ = value; }
  void setMaxIterations(int newMax)          { maxIters_ = newMax; }

  template<class SystemType>
  void solve(const SystemType & system, StateType & solutionInOut)
  {
    // deep copy the initial guess
    ::pressio::ops::deep_copy(this->template get<InitialGuessTag>(), solutionInOut);

    if (updateEnValue_ == Update::Standard)
    {
      auto extReg = reference_capture_registry_and_extend_with<
	StateTag, StateType &>(*this, solutionInOut);
      root_solving_loop_impl(tag_, system, extReg, stopEnValue_, stopTolerance_,
			     normDiagnostics_, diagnosticsLogger_, maxIters_,
			     DefaultUpdater());

    }
    else if (updateEnValue_ == Update::BacktrackStrictlyDecreasingObjective)
    {
      auto extReg = reference_capture_registry_and_extend_with<
	StateTag, LineSearchTrialStateTag,
	StateType &, StateType>(*this, solutionInOut, system.createState());

      root_solving_loop_impl(tag_, system, extReg, stopEnValue_, stopTolerance_,
			     normDiagnostics_, diagnosticsLogger_, maxIters_,
			     BacktrackStrictlyDecreasingObjectiveUpdater{});

    }
    else{
      throw std::runtime_error("Invalid criterion");
    }
  }

  // this method can be used when the solver is applied
  // to the same system used for constructing it
  void solve(StateType & solutionInOut)
  {
    auto * system = this->template get<SystemTag>();
    assert(system != nullptr);
    this->solve(*system, solutionInOut);
  }
};


}}}
#endif
