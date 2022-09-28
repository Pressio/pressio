
#ifndef PRESSIO_ROM_IMPL_GALERKIN_STEADY_MASKED_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_STEADY_MASKED_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
masked Galerkin problem represents:

   hypredOp masker[fom_r(phi x)] = 0

- fom_r is the fom "residual"
- phi is the basis
- masker is the masker

From this we get a "reduced" residual/jacobian:
R = phi^T hypredOp masker(fom_r(phi x))
J = phi^T hypredOp masker(dfom_r/dx(phi x) phi)
*/
template <
  class ReducedStateType,
  class ReducedResidualType,
  class ReducedJacobianType,
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HypRedOpType
  >
class GalerkinSteadyMaskedSystem
{

  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;

  // deduce the unmasked types
  using unmasked_fom_residual_type = typename FomSystemType::residual_type;
  using unmasked_fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<basis_matrix_type const &>()));

  // deduce the masked types
  using masked_fom_residual_type =
    decltype(std::declval<MaskerType const>().createApplyMaskResult
	     (std::declval<unmasked_fom_residual_type const &>()));

  using masked_fom_jac_action_result_type =
    decltype(std::declval<MaskerType const>().createApplyMaskResult
	     (std::declval<unmasked_fom_jac_action_result_type const &>()));

public:
  // required aliases
  using state_type    = ReducedStateType;
  using residual_type = ReducedResidualType;
  using jacobian_type = ReducedJacobianType;

  GalerkinSteadyMaskedSystem() = delete;
  GalerkinSteadyMaskedSystem(const TrialSubspaceType & trialSubspace,
			     const FomSystemType & fomSystem,
			     const MaskerType & masker,
			     const HypRedOpType & hyperReducer)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      hyperReducer_(hyperReducer),
      masker_(masker),
      unMaskedFomResidual_(fomSystem.createResidual()),
      unMaskedFomJacAction_(fomSystem.createApplyJacobianResult(trialSubspace_.get().basisOfTranslatedSpace())),
      maskedFomResidual_(masker.createApplyMaskResult(unMaskedFomResidual_)),
      maskedFomJacAction_(masker.createApplyMaskResult(unMaskedFomJacAction_))
  {}

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  residual_type createResidual() const{
    return impl::CreateGalerkinRhs<residual_type>()(trialSubspace_.get().dimension());
  }

  jacobian_type createJacobian() const{
    return impl::CreateGalerkinJacobian<jacobian_type>()(trialSubspace_.get().dimension());
  }

  void residualAndJacobian(const state_type & reducedState,
			   residual_type & reducedResidual,
			   jacobian_type & reducedJacobian,
			   bool computeJacobian) const
  {
    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    fomSystem_.get().residual(fomState_, unMaskedFomResidual_);
    masker_(unMaskedFomResidual_, maskedFomResidual_);
    hyperReducer_(maskedFomResidual_, reducedResidual);

    if (computeJacobian){
      fomSystem_.get().applyJacobian(fomState_, phi, unMaskedFomJacAction_);
      masker_(unMaskedFomJacAction_, maskedFomJacAction_);
      hyperReducer_(maskedFomJacAction_, reducedJacobian);
    }
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HypRedOpType> hyperReducer_;
  std::reference_wrapper<const MaskerType> masker_;

  // UNMASKED fom R,J instances
  mutable unmasked_fom_residual_type unMaskedFomResidual_;
  mutable unmasked_fom_jac_action_result_type unMaskedFomJacAction_;

  // MASKED fom R,J instances
  mutable masked_fom_residual_type maskedFomResidual_;
  mutable masked_fom_jac_action_result_type maskedFomJacAction_;
};

}}}
#endif
