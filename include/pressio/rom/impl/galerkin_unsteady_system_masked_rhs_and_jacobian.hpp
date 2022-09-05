
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_MASKED_RHS_AND_JACOBIAN_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_MASKED_RHS_AND_JACOBIAN_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  masked implicit galerkin system represents:

     d hat{y}/dt = hrOp masker ( fom_rhs(phi*hat{y}, ...) )

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
- hrOp is the hypred operator
so that it boils down to:

rhs = hrOp masked(fom_rhs(phi*hat{y}, ...))
rhs_jacobian = hrOp masked(d(fom_rhs(phi*hat{y}, ...))/dy phi)

*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class ReducedJacobianType,
  class TrialSpaceType,
  class FomSystemType,
  class RhsMaskerType,
  class JacobianActionMaskerType,
  class HyperReductionOperator
  >
class GalerkinMaskedOdeSystemRhsAndJacobian
{
  using basis_type = typename TrialSpaceType::basis_type;

  // deduce the unmasked types
  using unmasked_fom_rhs_type = typename FomSystemType::right_hand_side_type;
  using unmasked_fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<basis_type const &>()));

  // deduce the masked types
  using masked_fom_rhs_type = typename RhsMaskerType::result_type;
  using masked_fom_jac_action_result_type = typename JacobianActionMaskerType::result_type;

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using right_hand_side_type      = ReducedRhsType;
  using jacobian_type = ReducedJacobianType;

  GalerkinMaskedOdeSystemRhsAndJacobian() = delete;

  GalerkinMaskedOdeSystemRhsAndJacobian(const TrialSpaceType & trialSpace,
					const FomSystemType & fomSystem,
					const RhsMaskerType & rhsMasker,
					const JacobianActionMaskerType & jaMasker,
					const HyperReductionOperator & hrOp)
    : trialSpace_(trialSpace),
      fomSystem_(fomSystem),
      fomState_(trialSpace.createFullState()),
      hrOp_(hrOp),
      rhsMasker_(rhsMasker),
      jaMasker_(jaMasker),
      unMaskedFomRhs_(fomSystem.createRightHandSide()),
      unMaskedFomJacAction_(fomSystem.createApplyJacobianResult(trialSpace_.get().viewBasis())),
      maskedFomRhs_(rhsMasker.createApplyMaskResult(unMaskedFomRhs_)),
      maskedFomJacAction_(jaMasker.createApplyMaskResult(unMaskedFomJacAction_))
  {}

public:
  state_type createState() const{
    return trialSpace_.get().createReducedState();
  }

  right_hand_side_type createRightHandSide() const{
    const auto & phi = trialSpace_.get().viewBasis();
    return impl::CreateGalerkinRhs<right_hand_side_type>()(phi);
  }

  jacobian_type createJacobian() const{
    const auto & phi = trialSpace_.get().viewBasis();
    return impl::CreateGalerkinJacobian<jacobian_type>()(phi);
  }

  void rightHandSide(const state_type & reducedState,
		     const IndVarType & rhsEvaluationTime,
		     right_hand_side_type & reducedRhs) const
  {

    trialSpace_.get().mapFromReducedState(reducedState, fomState_);
    fomSystem_.get().rightHandSide(fomState_, rhsEvaluationTime, unMaskedFomRhs_);
    rhsMasker_(unMaskedFomRhs_, maskedFomRhs_);
    hrOp_(maskedFomRhs_, rhsEvaluationTime, reducedRhs);
  }

  void jacobian(const state_type & reducedState,
		const IndVarType & rhsEvaluationTime,
		jacobian_type & reducedJacobian) const

  {
    trialSpace_.get().mapFromReducedState(reducedState, fomState_);

    const auto & phi = trialSpace_.get().viewBasis();
    fomSystem_.get().applyJacobian(fomState_, phi, rhsEvaluationTime, unMaskedFomJacAction_);
    jaMasker_(unMaskedFomJacAction_, maskedFomJacAction_);
    hrOp_(maskedFomJacAction_, rhsEvaluationTime, reducedJacobian);
  }

protected:
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReductionOperator> hrOp_;

  // masker
  std::reference_wrapper<const RhsMaskerType> rhsMasker_;
  std::reference_wrapper<const JacobianActionMaskerType> jaMasker_;

  // UNMASKED objects
  mutable unmasked_fom_rhs_type unMaskedFomRhs_;
  mutable unmasked_fom_jac_action_result_type unMaskedFomJacAction_;

  // MASKED objects
  mutable masked_fom_rhs_type maskedFomRhs_;
  mutable masked_fom_jac_action_result_type maskedFomJacAction_;
};

}}} // end pressio::rom::impl
#endif
