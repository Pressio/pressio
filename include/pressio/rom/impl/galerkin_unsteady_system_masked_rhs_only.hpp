
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_MASKED_RHS_ONLY_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_MASKED_RHS_ONLY_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  masked explicit galerkin system represents:

     d hat{y}/dt = HrOp mask( fom_rhs(phi*hat{y}, ...) )

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
- HrOp is the hyper-red operator
*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class TrialSpaceType,
  class FomSystemType,
  class RhsMaskerType,
  class HyperReductionOperator
  >
class GalerkinMaskedOdeSystemOnlyRhs
{
  // deduce types
  using unmasked_fom_rhs_type = typename FomSystemType::right_hand_side_type;
  using masked_fom_rhs_type = typename RhsMaskerType::result_type;

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using right_hand_side_type      = ReducedRhsType;

  GalerkinMaskedOdeSystemOnlyRhs() = delete;

  GalerkinMaskedOdeSystemOnlyRhs(const TrialSpaceType & trialSpace,
				 const FomSystemType & fomSystem,
				 const RhsMaskerType & rhsMasker,
				 const HyperReductionOperator & hrOp)
    : trialSpace_(trialSpace),
      fomSystem_(fomSystem),
      fomState_(trialSpace.createFullState()),
      hrOp_(hrOp),
      rhsMasker_(rhsMasker),
      unMaskedFomRhs_(fomSystem.createRightHandSide()),
      maskedFomRhs_(rhsMasker.createApplyMaskResult(unMaskedFomRhs_))
  {}

public:
  state_type createState() const{
    return trialSpace_.get().createReducedState();
  }

  right_hand_side_type createRightHandSide() const{
    const auto & phi = trialSpace_.get().viewBasis();
    return impl::CreateGalerkinRhs<right_hand_side_type>()(phi);
  }

  void rightHandSide(const state_type & reducedState,
		     const IndVarType & rhsEvaluationTime,
		     right_hand_side_type & reducedRhs) const
  {

    // reconstruct fom state fomState = phi*reducedState
    trialSpace_.get().mapFromReducedState(reducedState, fomState_);

    // evaluate fomRhs and mask it
    fomSystem_.get().rightHandSide(fomState_, rhsEvaluationTime, unMaskedFomRhs_);
    rhsMasker_(unMaskedFomRhs_, maskedFomRhs_);

    // evaluate reduced rhs
    hrOp_(maskedFomRhs_, rhsEvaluationTime, reducedRhs);
  }

protected:
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReductionOperator> hrOp_;

  std::reference_wrapper<const RhsMaskerType> rhsMasker_;
  mutable unmasked_fom_rhs_type unMaskedFomRhs_;
  mutable masked_fom_rhs_type maskedFomRhs_;
};

}}} // end pressio::rom::impl
#endif
