
#ifndef ROM_IMPL_GALERKIN_UNSTEADY_EXPLICIT_PROBLEM_HPP_
#define ROM_IMPL_GALERKIN_UNSTEADY_EXPLICIT_PROBLEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <class GalSystem>
class GalerkinUnsteadyExplicitProblem
{
  // note: to deduce the stepper_type it does not really matter
  // what scheme enum value we use, as long as it is an explicit one
  using stepper_type =
    decltype(::pressio::ode::create_explicit_stepper
	     (::pressio::ode::StepScheme::ForwardEuler,
	      std::declval<GalSystem &>()
	      ));

public:
  // required aliases to be steppable
  using state_type = typename GalSystem::state_type;
  using independent_variable_type  = typename GalSystem::independent_variable_type;

  template<class ...Args>
  GalerkinUnsteadyExplicitProblem(::pressio::ode::StepScheme schemeName,
				  Args && ... args)
    : galSystem_(std::forward<Args>(args)...),
      stepper_( ::pressio::ode::create_explicit_stepper(schemeName, galSystem_) )
  {}

  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> sStart,
		  pressio::ode::StepCount sCount,
		  pressio::ode::StepSize<independent_variable_type> sSize)
  {
    stepper_(state, sStart, sCount, sSize);
  }

private:
  GalSystem galSystem_;
  stepper_type stepper_;
};


/*
  default explicit galerkin system with time and state dependent
  mass matrix represents:

     [phi^T MM(y, time, ...) phi] d hat{y}/dt = phi^T fom_rhs(phi*hat{y}, ...)

- y is the fom state
- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
- M_gal = [phi^T MM(y, time, ...) phi] is the reduced mass matrix
*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedMassMatrixType,
  class TrialSubspaceType,
  class FomMassMatrixOperatorType
  >
class GalerkinExplicitReducedVariableMassMatrixEvaluator
{

  // deduce from the fom object the type of result of
  // applying the mass matrix to the basis
  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;
  using fom_mass_matrix_action_result_type =
    decltype(std::declval<FomMassMatrixOperatorType const>().createApplyMassMatrixResult
	     (std::declval<basis_matrix_type const &>()) );

public:
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using mass_matrix_type          = ReducedMassMatrixType;

  GalerkinExplicitReducedVariableMassMatrixEvaluator(const TrialSubspaceType & trialSubspace,
					     const FomMassMatrixOperatorType & fomMassMatOp)
    : trialSubspace_(trialSubspace),
      fomState_(trialSubspace.createFullState()),
      fomMassMatOp_(fomMassMatOp),
      fomMassMatAction_(fomMassMatOp.createApplyMassMatrixResult(trialSubspace_.get().basisOfTranslatedSpace()))
  {}

  mass_matrix_type createMassMatrix() const{
    return impl::CreateGalerkinMassMatrix<mass_matrix_type>()(trialSubspace_.get().dimension());
  }

  void massMatrix(const state_type & reducedState,
		  const independent_variable_type & evalTime,
		  mass_matrix_type & reducedMassMatrix) const
  {
    // reconstruct fom state fomState = phi*reducedState
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();

    // evaluate fom mass matrix action
    fomMassMatOp_(fomState_, phi, evalTime, fomMassMatAction_);

    // compute the reduced mass matrix
    using phi_scalar_t   = typename ::pressio::Traits<basis_matrix_type>::scalar_type;
    constexpr auto alpha = ::pressio::utils::Constants<phi_scalar_t>::one();
    using rmm_scalar_t  = typename ::pressio::Traits<mass_matrix_type>::scalar_type;
    constexpr auto beta = ::pressio::utils::Constants<rmm_scalar_t>::zero();
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
			    alpha, phi, fomMassMatAction_,
			    beta, reducedMassMatrix);
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  mutable typename TrialSubspaceType::full_state_type fomState_;
  std::reference_wrapper<const FomMassMatrixOperatorType> fomMassMatOp_;
  mutable fom_mass_matrix_action_result_type fomMassMatAction_;
};

template <
  class GalSystem,
  class ReducedMassMatrixOperator
  >
class GalerkinUnsteadyWithMassMatrixExplicitProblem
{
  // note: to deduce the stepper_type it does not really matter
  // what scheme enum value we use, as long as it is an explicit one
  using stepper_type =
    decltype(::pressio::ode::create_explicit_stepper
	     (::pressio::ode::StepScheme::ForwardEuler,
	      std::declval<GalSystem &>(),
	      std::declval<ReducedMassMatrixOperator &>()
	      ));

public:
  // required aliases to be steppable
  using state_type = typename GalSystem::state_type;
  using independent_variable_type = typename GalSystem::independent_variable_type;
  using mass_matrix_type = typename ReducedMassMatrixOperator::mass_matrix_type;

  template<class TrialSubspaceType, class FomSystemType, class FomMassMatrixOperatorType>
  GalerkinUnsteadyWithMassMatrixExplicitProblem(::pressio::ode::StepScheme schemeName,
						const TrialSubspaceType & trialSubspace,
						const FomSystemType & fomSystem,
						const FomMassMatrixOperatorType & fomMMOp)
    : galSystem_(trialSubspace, fomSystem),
      galMassMat_(trialSubspace, fomMMOp),
      stepper_( ::pressio::ode::create_explicit_stepper(schemeName, galSystem_, galMassMat_) )
  {}

  template<class LinearSolverType>
  void operator()(state_type & state,
		  pressio::ode::StepStartAt<independent_variable_type> sStart,
		  pressio::ode::StepCount sCount,
		  pressio::ode::StepSize<independent_variable_type> sSize,
		  LinearSolverType & linSolver)
  {
    stepper_(state, sStart, sCount, sSize, linSolver);
  }

private:
  GalSystem galSystem_;
  ReducedMassMatrixOperator galMassMat_;
  stepper_type stepper_;
};

}}} // end pressio::rom::impl
#endif  // ROM_IMPL_GALERKIN_UNSTEADY_EXPLICIT_PROBLEM_HPP_
