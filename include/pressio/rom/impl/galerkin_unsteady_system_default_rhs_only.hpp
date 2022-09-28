
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_DEFAULT_RHS_ONLY_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_DEFAULT_RHS_ONLY_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  default explicit galerkin system represents:

     d hat{y}/dt = phi^T fom_rhs(phi*hat{y}, ...)

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis

*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class TrialSubspaceType,
  class FomSystemType
  >
class GalerkinDefaultOdeSystemOnlyRhs
{
  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using right_hand_side_type      = ReducedRhsType;

  GalerkinDefaultOdeSystemOnlyRhs() = delete;

  GalerkinDefaultOdeSystemOnlyRhs(const TrialSubspaceType & trialSubspace,
				  const FomSystemType & fomSystem)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      fomRhs_(fomSystem.createRightHandSide())
  {}

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  right_hand_side_type createRightHandSide() const{
    return impl::CreateGalerkinRhs<right_hand_side_type>()(trialSubspace_.get().dimension());
  }

  void rightHandSide(const state_type & reducedState,
		     const IndVarType & rhsEvaluationTime,
		     right_hand_side_type & reducedRhs) const
  {

    // reconstruct fom state fomState = phi*reducedState
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    // evaluate fomRhs
    fomSystem_.get().rightHandSide(fomState_, rhsEvaluationTime, fomRhs_);

    // compute the reduced rhs
    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    using phi_scalar_t = typename ::pressio::Traits<basis_matrix_type>::scalar_type;
    constexpr auto alpha = ::pressio::utils::Constants<phi_scalar_t>::one();
    using rhs_scalar_t = typename ::pressio::Traits<right_hand_side_type>::scalar_type;
    constexpr auto beta = ::pressio::utils::Constants<rhs_scalar_t>::zero();
    ::pressio::ops::product(::pressio::transpose(),
			    alpha, phi, fomRhs_,
			    beta, reducedRhs);
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  mutable typename FomSystemType::right_hand_side_type fomRhs_;
};

}}} // end pressio::rom::impl
#endif
