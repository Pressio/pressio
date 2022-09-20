
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_HYPRED_RHS_ONLY_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_HYPRED_RHS_ONLY_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  hyp-red explicit galerkin system represents:

     d hat{y}/dt = HrOp fom_rhs(phi*hat{y}, ...)

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
  class HyperReductionOperator
  >
class GalerkinHyperReducedOdeSystemOnlyRhs
{

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using right_hand_side_type      = ReducedRhsType;

  GalerkinHyperReducedOdeSystemOnlyRhs() = delete;

  GalerkinHyperReducedOdeSystemOnlyRhs(const TrialSpaceType & trialSpace,
				       const FomSystemType & fomSystem,
				       const HyperReductionOperator & hrOp)
    : trialSpace_(trialSpace),
      fomSystem_(fomSystem),
      fomState_(trialSpace.createFullState()),
      hrOp_(hrOp),
      fomRhs_(fomSystem.createRightHandSide())
  {}

public:
  state_type createState() const{
    return trialSpace_.get().createReducedState();
  }

  right_hand_side_type createRightHandSide() const{
    const auto & phi = trialSpace_.get().basisOfTranslatedSpace();
    return impl::CreateGalerkinRhs<right_hand_side_type>()(phi);
  }

  void rightHandSide(const state_type & reducedState,
		     const IndVarType & rhsEvaluationTime,
		     right_hand_side_type & reducedRhs) const
  {

    // reconstruct fom state fomState = phi*reducedState
    trialSpace_.get().mapFromReducedState(reducedState, fomState_);

    // evaluate fomRhs
    fomSystem_.get().rightHandSide(fomState_, rhsEvaluationTime, fomRhs_);

    // evaluate reduced rhs
    hrOp_(fomRhs_, rhsEvaluationTime, reducedRhs);
  }

private:
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReductionOperator> hrOp_;
  mutable typename FomSystemType::right_hand_side_type fomRhs_;
};

}}} // end pressio::rom::impl
#endif
