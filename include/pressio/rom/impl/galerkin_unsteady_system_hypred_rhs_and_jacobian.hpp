
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_HYPRED_RHS_AND_JAC_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_HYPRED_RHS_AND_JAC_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  hypred implicit galerkin system represents:

     d hat{y}/dt = hrOp fom_rhs(phi*hat{y}, ...)

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
- hrOp is the hypred operator
so that it boils down to:

rhs = hrOp fom_rhs(phi*hat{y}, ...)
rhs_jacobian = hrOp d(fom_rhs(phi*hat{y}, ...))/dy phi

*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class ReducedJacobianType,
  class TrialSpaceType,
  class FomSystemType,
  class HyperReductionOperator
  >
class GalerkinHypRedOdeSystemRhsAndJacobian
{

  // deduce from the fom object the type of result of
  // applying the Jacobian to the basis
  using basis_matrix_type = typename TrialSpaceType::basis_matrix_type;
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

  GalerkinHypRedOdeSystemRhsAndJacobian(const TrialSpaceType & trialSpace,
					const FomSystemType & fomSystem,
					const HyperReductionOperator & hrOp)
    : trialSpace_(trialSpace),
      fomSystem_(fomSystem),
      fomState_(trialSpace.createFullState()),
      hrOp_(hrOp),
      fomRhs_(fomSystem.createRightHandSide()),
      fomJacAction_(fomSystem.createApplyJacobianResult(trialSpace_.get().basisOfTranslatedSpace()))
  {}

public:
  state_type createState() const{
    return trialSpace_.get().createReducedState();
  }

  right_hand_side_type createRightHandSide() const{
    const auto & phi = trialSpace_.get().basisOfTranslatedSpace();
    return impl::CreateGalerkinRhs<right_hand_side_type>()(phi);
  }

  jacobian_type createJacobian() const{
    const auto & phi = trialSpace_.get().basisOfTranslatedSpace();
    return impl::CreateGalerkinJacobian<jacobian_type>()(phi);
  }

  void rightHandSide(const state_type & reducedState,
		     const IndVarType & rhsEvaluationTime,
		     right_hand_side_type & reducedRhs) const
  {

    // reconstruct fom state fomState = phi*reducedState
    trialSpace_.get().mapFromReducedState(reducedState, fomState_);

    // evaluate fomRhs and apply hrOp
    fomSystem_.get().rightHandSide(fomState_, rhsEvaluationTime, fomRhs_);

    // compute the reduced rhs
    hrOp_(fomRhs_, rhsEvaluationTime, reducedRhs);
  }

  void jacobian(const state_type & reducedState,
		const IndVarType & rhsEvaluationTime,
		jacobian_type & reducedJacobian) const

  {
    // reconstruct fom state fomState = phi*reducedState
    trialSpace_.get().mapFromReducedState(reducedState, fomState_);

    const auto & phi = trialSpace_.get().basisOfTranslatedSpace();

    // evaluate fom jacobian action: fomJacAction_ = fom_J * phi
    fomSystem_.get().applyJacobian(fomState_, phi, rhsEvaluationTime, fomJacAction_);

    // compute the reduced jacobian
    hrOp_(fomJacAction_, rhsEvaluationTime, reducedJacobian);
  }

private:
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReductionOperator> hrOp_;
  mutable typename FomSystemType::right_hand_side_type fomRhs_;
  mutable fom_jac_action_result_type fomJacAction_;
};

}}} // end pressio::rom::impl
#endif
