
#ifndef PRESSIO_ROM_IMPL_LSPG_STEADY_MASKED_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_LSPG_STEADY_MASKED_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
LSPG steady masked represents:

  min_x || mask[fom_r(phi x)]||

- fom_r is the fom "residual"
- phi is the basis
*/
template <
  class ReducedStateType,
  class TrialSpaceType,
  class FomSystemType,
  class ResidualMaskerType,
  class JacobianActionMaskerType
  >
class LspgSteadyMaskedSystem
{

  using unmasked_fom_residual_type = typename FomSystemType::residual_type;
  using unmasked_fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
      (std::declval<typename TrialSpaceType::basis_matrix_type const &>())
      );

  // deduce the masked types
  using masked_fom_residual_type = typename ResidualMaskerType::result_type;
  using masked_fom_jac_action_result_type = typename JacobianActionMaskerType::result_type;

public:
  // required aliases
  using state_type    = ReducedStateType;
  using residual_type = masked_fom_residual_type;
  using jacobian_type = masked_fom_jac_action_result_type;

  LspgSteadyMaskedSystem() = delete;

  LspgSteadyMaskedSystem(TrialSpaceType & trialSpace,
			 const FomSystemType & fomSystem,
			 const ResidualMaskerType & rMasker,
			 const JacobianActionMaskerType & jaMasker)
    : trialSpace_(trialSpace),
      fomSystem_(fomSystem),
      fomState_(trialSpace.createFullState()),
      rMasker_(rMasker),
      jaMasker_(jaMasker),
      unMaskedFomResidual_(fomSystem.createResidual()),
      unMaskedFomJacAction_(fomSystem.createApplyJacobianResult(trialSpace_.get().basisOfTranslatedSpace()))
  {}

public:
  state_type createState() const{
    return trialSpace_.get().createReducedState();
  }

  residual_type createResidual() const{
    auto tmp = fomSystem_.get().createResidual();
    return rMasker_.get().createApplyMaskResult(tmp);
  }

  jacobian_type createJacobian() const{
    auto tmp = fomSystem_.get().createApplyJacobianResult(trialSpace_.get().basisOfTranslatedSpace());
    return jaMasker_.get().createApplyMaskResult(tmp);
  }

  void residualAndJacobian(const state_type & reducedState,
			   residual_type & lsgpResidual,
			   jacobian_type & lspgJacobian,
			   bool computeJacobian) const
  {
    trialSpace_.get().mapFromReducedState(reducedState, fomState_);

    fomSystem_.get().residual(fomState_, unMaskedFomResidual_);
    rMasker_(unMaskedFomResidual_, lsgpResidual);

    if (computeJacobian){
      const auto & phi = trialSpace_.get().basisOfTranslatedSpace();
      fomSystem_.get().applyJacobian(fomState_, phi, unMaskedFomJacAction_);
      jaMasker_(unMaskedFomJacAction_, lspgJacobian);
    }
  }

protected:
  std::reference_wrapper<TrialSpaceType> trialSpace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;

  // masker
  std::reference_wrapper<const ResidualMaskerType> rMasker_;
  std::reference_wrapper<const JacobianActionMaskerType> jaMasker_;

  // UNMASKED fom R,J instances
  mutable unmasked_fom_residual_type unMaskedFomResidual_;
  mutable unmasked_fom_jac_action_result_type unMaskedFomJacAction_;
};

}}} // end pressio::rom::impl
#endif
