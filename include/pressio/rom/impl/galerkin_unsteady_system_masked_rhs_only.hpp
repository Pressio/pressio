
#ifndef ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_MASKED_RHS_ONLY_HPP_
#define ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_MASKED_RHS_ONLY_HPP_

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
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType
  >
class GalerkinMaskedOdeSystemOnlyRhs
{
  // deduce types
  using unmasked_fom_rhs_type = typename FomSystemType::right_hand_side_type;
  using masked_fom_rhs_type =
    decltype(std::declval<MaskerType const>().createResultOfMaskActionOn
	     (std::declval<unmasked_fom_rhs_type const &>()));

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using right_hand_side_type      = ReducedRhsType;

  GalerkinMaskedOdeSystemOnlyRhs() = delete;

  GalerkinMaskedOdeSystemOnlyRhs(const TrialSubspaceType & trialSubspace,
				 const FomSystemType & fomSystem,
				 const MaskerType & masker,
				 const HyperReducerType & hyperReducer)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      hyperReducer_(hyperReducer),
      masker_(masker),
      unMaskedFomRhs_(fomSystem.createRightHandSide()),
      maskedFomRhs_(masker.createResultOfMaskActionOn(unMaskedFomRhs_))
  {}

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  right_hand_side_type createRightHandSide() const{
    return impl::CreateGalerkinRhs<right_hand_side_type>()(trialSubspace_.get().dimension());
  }

  void operator()(const state_type & reducedState,
		  const IndVarType & rhsEvaluationTime,
		  right_hand_side_type & reducedRhs) const
  {

    // reconstruct fom state fomState = phi*reducedState
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    // evaluate fomRhs and mask it
    fomSystem_.get().rightHandSide(fomState_, rhsEvaluationTime, unMaskedFomRhs_);
    masker_(unMaskedFomRhs_, maskedFomRhs_);

    // evaluate reduced rhs
    hyperReducer_(maskedFomRhs_, rhsEvaluationTime, reducedRhs);
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReducerType> hyperReducer_;
  std::reference_wrapper<const MaskerType> masker_;
  mutable unmasked_fom_rhs_type unMaskedFomRhs_;
  mutable masked_fom_rhs_type maskedFomRhs_;
};

}}} // end pressio::rom::impl
#endif  // ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_MASKED_RHS_ONLY_HPP_
