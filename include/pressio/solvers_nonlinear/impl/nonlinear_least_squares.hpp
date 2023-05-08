
#ifndef SOLVERS_NONLINEAR_IMPL_LEAST_SQUARES_SOLVER_HPP_
#define SOLVERS_NONLINEAR_IMPL_LEAST_SQUARES_SOLVER_HPP_

namespace pressio{
namespace nonlinearsolvers{
namespace impl{

template<
  class ProblemTag,
  class UserDefinedSystemType,
  class RegistryType,
  class ToleranceType,
  class DiagnosticsContainerType,
  class DiagnosticsLoggerType,
  class UpdaterType>
void nonlin_ls_solving_loop_impl(ProblemTag problemTag,
				 const UserDefinedSystemType & system,
				 RegistryType & reg,
				 Stop stopEnumValue,
				 ToleranceType stopTolerance,
				 DiagnosticsContainerType & normDiagnostics,
				 const DiagnosticsLoggerType & logger,
				 int maxIters,
				 UpdaterType && updater)
{

  auto mustStop = [
#ifdef PRESSIO_ENABLE_CXX17
		   &normDiag = std::as_const(normDiagnostics),
#else
		   &normDiag = normDiagnostics,
#endif
		   stopEnumValue, maxIters, stopTolerance](int stepCount)
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
    const bool isFirstIteration = iStep==1;

    // 1. compute operators
    try{
      auto objValue = compute_nonlinearls_operators_and_objective(problemTag, reg, system);
      normDiagnostics[InternalDiagnostic::objectiveAbsoluteRelative].update(objValue, isFirstIteration);
    }
    catch (::pressio::eh::ResidualEvaluationFailureUnrecoverable const &e){
      PRESSIOLOG_CRITICAL(e.what());
      throw ::pressio::eh::NonlinearSolveFailure();
    }
    catch (::pressio::eh::ResidualHasNans const &e){
      PRESSIOLOG_CRITICAL(e.what());
      throw ::pressio::eh::NonlinearSolveFailure();
    }

    // 2. solve for correction
    compute_correction(problemTag, reg);

    /* stage 3 */
    std::for_each(normDiagnostics.begin(), normDiagnostics.end(),
		  [&reg, isFirstIteration](auto & v){
		    compute_norm_internal_diagnostics(reg, isFirstIteration, v);
		  });
    logger(iStep, normDiagnostics);

    /* stage 4*/
    if (mustStop(iStep)){
      PRESSIOLOG_DEBUG("nonlinsolver: stopping");
      break;
    }

    // 5. run update and continue
    try{
      using state_type = typename UserDefinedSystemType::state_type;
      auto objective = [&reg, &system, problemTag](const state_type & stateIn){
	return compute_nonlinearls_objective(problemTag, reg, stateIn, system);
      };
      const auto currObjValue =
	normDiagnostics[InternalDiagnostic::objectiveAbsoluteRelative].getAbsolute();
      updater(reg, objective, currObjValue);

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


template<class Tag, class StateType, class RegistryType, class ScalarType>
class NonLinLeastSquares
{
  static_assert(std::is_floating_point<ScalarType>::value,
		"Impl currently only supporting floating point");

  Tag tag_;
  int maxIters_ = 100;
  Stop stopEnValue_ = Stop::WhenAbsolutel2NormOfCorrectionBelowTolerance;
  ScalarType stopTolerance_ = static_cast<ScalarType>(0.000001);
  Update updateEnValue_ = Update::Standard;

  RegistryType reg_;

  using diagnostics_container = DiagnosticsContainer<
    InternalDiagnosticDataWithAbsoluteRelativeTracking<ScalarType> >;
  diagnostics_container diagnostics_;
  DiagnosticsLogger diagnosticsLogger_ = {};

public:
  NonLinLeastSquares(Tag tagIn,
		     RegistryType && reg,
		     const std::vector<Diagnostic> & diags)
    : tag_(tagIn), reg_(std::move(reg)), diagnostics_(diags)
  {

    // currently we don't have the diagonostics stuff all flushed out
    // so we limit it to work for a specific case
    const auto & publicDiags = diagnostics_.publicNames();
    assert(publicDiags.size() == 6);
    assert(publicDiags[0] == Diagnostic::objectiveAbsolute);
    assert(publicDiags[1] == Diagnostic::objectiveRelative);
    assert(publicDiags[2] == Diagnostic::correctionAbsolutel2Norm);
    assert(publicDiags[3] == Diagnostic::correctionRelativel2Norm);
    assert(publicDiags[4] == Diagnostic::gradientAbsolutel2Norm);
    assert(publicDiags[5] == Diagnostic::gradientRelativel2Norm);

    // check that the stopping criterion uses a metric already
    // supported in the diagonstics, otherwise add it
    const auto stopMetric = stop_criterion_to_public_diagnostic(stopEnValue_);
    diagnostics_.addIfUnsupported(stopMetric);
    // need to reset the logger since the names might have changed
    diagnosticsLogger_.resetFor(publicDiags);
  }

  // query/set update criterion
  Update currentUpdateCriterion() const  { return updateEnValue_; }
  void setUpdateCriterion(Update value)  { updateEnValue_ = value; }

  // query/set stop criterion, tolerance
  Stop currentStopCriterion() const       { return stopEnValue_; }
  void setStopCriterion(Stop value)	  { stopEnValue_ = value; }
  void setStopTolerance(ScalarType value) { stopTolerance_ = value; }
  void setMaxIterations(int newMax)       { maxIters_ = newMax; }

  // this method can be used when passing a system object
  // that is different but syntactically and semantically equivalent
  // to the one used for constructing the solver
  template<class SystemType>
  void solve(const SystemType & system, StateType & solutionInOut)
  {
    // the solve method potentially be called multiple times
    // so we need to reset the data in the registry everytime
    // if not applicable, this is a noop
    reset_for_new_solve_loop(tag_, reg_);

    switch (updateEnValue_)
    {
      case Update::Standard:
	solve_with_standard_update_impl(system, solutionInOut);
	break;

      case Update::BacktrackStrictlyDecreasingObjective:
      case Update::Armijo:{
	solve_with_line_search_impl(system, solutionInOut);
	break;
      }

      case Update::LMSchedule1:
      case Update::LMSchedule2:{
	solve_lm_impl(system, solutionInOut);
	break;
      }

      default:
	throw std::runtime_error("invalid Update value");
    };
  }

  // this method can be used when the solver is applied
  // to the same system used for constructing it
  void solve(StateType & solutionInOut)
  {
    auto * system = reg_.template get<SystemTag>();
    assert(system != nullptr);
    this->solve(*system, solutionInOut);
  }

private:
  template<class SystemType>
  void solve_with_standard_update_impl(const SystemType & system,
				       StateType & solutionInOut)
  {
    auto extReg = reference_capture_registry_and_extend_with<
      StateTag, StateType &>(reg_, solutionInOut);

    // the solve method potentially be called multiple times
    // so we need to reset the data in the registry everytime
    reset_for_new_solve_loop(tag_, extReg);

    nonlin_ls_solving_loop_impl(tag_, system, extReg,
				stopEnValue_, stopTolerance_,
				diagnostics_, diagnosticsLogger_,
				maxIters_,
				DefaultUpdater());
  }

  template<class SystemType>
  void solve_with_line_search_impl(const SystemType & system,
				   StateType & solutionInOut)
  {
    auto extReg = reference_capture_registry_and_extend_with<
      StateTag, LineSearchTrialStateTag,
      StateType &, StateType>(reg_, solutionInOut, system.createState());

    // the solve method potentially be called multiple times
    // so we need to reset the data in the registry everytime
    reset_for_new_solve_loop(tag_, extReg);

    if (updateEnValue_ == Update::BacktrackStrictlyDecreasingObjective){
      nonlin_ls_solving_loop_impl(tag_, system, extReg,
				  stopEnValue_, stopTolerance_,
				  diagnostics_, diagnosticsLogger_,
				  maxIters_,
				  BacktrackStrictlyDecreasingObjectiveUpdater{});
    }else{
      nonlin_ls_solving_loop_impl(tag_, system, extReg,
				  stopEnValue_, stopTolerance_,
				  diagnostics_, diagnosticsLogger_,
				  maxIters_,
				  BacktrackStrictlyDecreasingObjectiveUpdater{});
    }
  }

  template<class SystemType, class _Tag = Tag>
  mpl::enable_if_t< std::is_same<_Tag, LevenbergMarquardtNormalEqTag>::value >
  solve_lm_impl(const SystemType & system,
		StateType & solutionInOut)
  {
    auto extReg = reference_capture_registry_and_extend_with<
      StateTag, LineSearchTrialStateTag,
      StateType &, StateType>(reg_, solutionInOut, system.createState());

    if (updateEnValue_ == Update::LMSchedule1){
      using up_t = LMSchedule1Updater<ScalarType, StateType>;
      nonlin_ls_solving_loop_impl(tag_, system, extReg,
				  stopEnValue_, stopTolerance_,
				  diagnostics_, diagnosticsLogger_,
				  maxIters_,
				  up_t(solutionInOut));
    }else{
      using up_t = LMSchedule2Updater<ScalarType, StateType>;
      nonlin_ls_solving_loop_impl(tag_, system, extReg,
				  stopEnValue_, stopTolerance_,
				  diagnostics_, diagnosticsLogger_,
				  maxIters_,
				  up_t(solutionInOut));
    }
  }

  template<class SystemType, class _Tag = Tag>
  mpl::enable_if_t< !std::is_same<_Tag, LevenbergMarquardtNormalEqTag>::value >
  solve_lm_impl(const SystemType & /*system*/,
		StateType & /*solutionInOut*/)
  {
    throw std::runtime_error("invalid Update value");
  }
};

}}}
#endif
