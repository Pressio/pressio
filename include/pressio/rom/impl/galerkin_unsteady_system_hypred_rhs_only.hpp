
#ifndef ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_HYPRED_RHS_ONLY_HPP_
#define ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_HYPRED_RHS_ONLY_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  hyp-red explicit galerkin system represents:

     d hat{y}/dt = hyperReducer fom_rhs(phi*hat{y}, ...)

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
- hyperReducer is the hyper-red operator
*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType
  >
class GalerkinHyperReducedOdeSystemOnlyRhs
{

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using right_hand_side_type      = ReducedRhsType;

  GalerkinHyperReducedOdeSystemOnlyRhs() = delete;

  GalerkinHyperReducedOdeSystemOnlyRhs(const TrialSubspaceType & trialSubspace,
				       const FomSystemType & fomSystem,
				       const HyperReducerType & hyperReducer)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      hyperReducer_(hyperReducer),
      fomRhs_(fomSystem.createRightHandSide())
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

    // evaluate fomRhs
    fomSystem_.get().rightHandSide(fomState_, rhsEvaluationTime, fomRhs_);

    // evaluate reduced rhs
    hyperReducer_(fomRhs_, rhsEvaluationTime, reducedRhs);
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReducerType> hyperReducer_;
  mutable typename FomSystemType::right_hand_side_type fomRhs_;
};

}}} // end pressio::rom::impl
#endif  // ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_HYPRED_RHS_ONLY_HPP_
