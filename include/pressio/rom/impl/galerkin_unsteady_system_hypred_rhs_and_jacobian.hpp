
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_HYPRED_RHS_AND_JAC_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_HYPRED_RHS_AND_JAC_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  hypred implicit galerkin system represents:

     d hat{y}/dt = hyperReducer fom_rhs(phi*hat{y}, ...)

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
- hyperReducer is the hypred operator
so that it boils down to:

rhs = hyperReducer fom_rhs(phi*hat{y}, ...)
rhs_jacobian = hyperReducer d(fom_rhs(phi*hat{y}, ...))/dy phi

*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class ReducedJacobianType,
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType
  >
class GalerkinHypRedOdeSystemRhsAndJacobian
{

  // deduce from the fom object the type of result of
  // applying the Jacobian to the basis
  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<basis_matrix_type const &>()) );

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using right_hand_side_type      = ReducedRhsType;
  using jacobian_type = ReducedJacobianType;

  GalerkinHypRedOdeSystemRhsAndJacobian() = delete;

  GalerkinHypRedOdeSystemRhsAndJacobian(const TrialSubspaceType & trialSubspace,
					const FomSystemType & fomSystem,
					const HyperReducerType & hyperReducer)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      hyperReducer_(hyperReducer),
      fomRhs_(fomSystem.createRightHandSide()),
      fomJacAction_(fomSystem.createApplyJacobianResult(trialSubspace_.get().basisOfTranslatedSpace()))
  {}

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  right_hand_side_type createRightHandSide() const{
    return impl::CreateGalerkinRhs<right_hand_side_type>()(trialSubspace_.get().dimension());
  }

  jacobian_type createJacobian() const{
    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    return impl::CreateGalerkinJacobian<jacobian_type>()(trialSubspace_.get().dimension());
  }

  void rightHandSide(const state_type & reducedState,
		     const IndVarType & rhsEvaluationTime,
		     right_hand_side_type & reducedRhs) const
  {

    // reconstruct fom state fomState = phi*reducedState
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    // evaluate fomRhs and apply hyperReducer
    fomSystem_.get().rightHandSide(fomState_, rhsEvaluationTime, fomRhs_);

    // compute the reduced rhs
    hyperReducer_(fomRhs_, rhsEvaluationTime, reducedRhs);
  }

  void jacobian(const state_type & reducedState,
		const IndVarType & rhsEvaluationTime,
		jacobian_type & reducedJacobian) const

  {
    // reconstruct fom state fomState = phi*reducedState
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();

    // evaluate fom jacobian action: fomJacAction_ = fom_J * phi
    fomSystem_.get().applyJacobian(fomState_, phi, rhsEvaluationTime, fomJacAction_);

    // compute the reduced jacobian
    hyperReducer_(fomJacAction_, rhsEvaluationTime, reducedJacobian);
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReducerType> hyperReducer_;
  mutable typename FomSystemType::right_hand_side_type fomRhs_;
  mutable fom_jac_action_result_type fomJacAction_;
};

}}} // end pressio::rom::impl
#endif
