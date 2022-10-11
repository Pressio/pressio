
#ifndef ROM_IMPL_LSPG_STEADY_SYSTEM_MASKED_HPP_
#define ROM_IMPL_LSPG_STEADY_SYSTEM_MASKED_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
LSPG steady masked represents:

  min_x || mask[fom_r(phi x)]||

- fom_r is the fom "residual"
- phi is the basis
*/
template <
  class ReducedStateType,
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType
  >
class LspgSteadyMaskedSystem
{

  using unmasked_fom_residual_type = typename FomSystemType::residual_type;
  using unmasked_fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createResultOfJacobianActionOn
      (std::declval<typename TrialSubspaceType::basis_matrix_type const &>())
      );

  // deduce the masked types
  using masked_fom_residual_type =
    decltype(std::declval<MaskerType const>().createResultOfMaskActionOn
	     (std::declval<unmasked_fom_residual_type const &>()));

  using masked_fom_jac_action_result_type =
    decltype(std::declval<MaskerType const>().createResultOfMaskActionOn
	     (std::declval<unmasked_fom_jac_action_result_type const &>()));

public:
  // required aliases
  using state_type    = ReducedStateType;
  using residual_type = masked_fom_residual_type;
  using jacobian_type = masked_fom_jac_action_result_type;

  LspgSteadyMaskedSystem(TrialSubspaceType & trialSubspace,
			 const FomSystemType & fomSystem,
			 const MaskerType & masker)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      masker_(masker),
      unMaskedFomResidual_(fomSystem.createResidual()),
      unMaskedFomJacAction_(fomSystem.createResultOfJacobianActionOn(trialSubspace_.get().basisOfTranslatedSpace()))
  {}

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  residual_type createResidual() const{
    auto tmp = fomSystem_.get().createResidual();
    return masker_.get().createResultOfMaskActionOn(tmp);
  }

  jacobian_type createJacobian() const{
    auto tmp = fomSystem_.get().createResultOfJacobianActionOn(trialSubspace_.get().basisOfTranslatedSpace());
    return masker_.get().createResultOfMaskActionOn(tmp);
  }

  void residualAndJacobian(const state_type & reducedState,
			   residual_type & lsgpResidual,
			   jacobian_type & lspgJacobian,
			   bool computeJacobian) const
  {
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    fomSystem_.get().residualAndJacobianAction(fomState_,
					       unMaskedFomResidual_,
					       phi, unMaskedFomJacAction_,
					       computeJacobian);
    masker_(unMaskedFomResidual_, lsgpResidual);
    if (computeJacobian){
      masker_(unMaskedFomJacAction_, lspgJacobian);
    }
  }

protected:
  std::reference_wrapper<TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const MaskerType> masker_;
  mutable unmasked_fom_residual_type unMaskedFomResidual_;
  mutable unmasked_fom_jac_action_result_type unMaskedFomJacAction_;
};

}}} // end pressio::rom::impl
#endif  // ROM_IMPL_LSPG_STEADY_SYSTEM_MASKED_HPP_
