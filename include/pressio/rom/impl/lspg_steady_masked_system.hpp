
#ifndef PRESSIO_ROM_IMPL_LSPG_STEADY_MASKED_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_LSPG_STEADY_MASKED_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  class LspgStateType,
  class TrialSpaceType,
  class FomSystemType,
  class MaskerType
  >
class LspgSteadyMaskedSystem
{

  using unmasked_residual_type = typename FomSystemType::residual_type;
  using unmasked_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
      (std::declval<typename TrialSpaceType::basis_type const &>())
      );

  using masked_residual_type =
    decltype(std::declval<MaskerType const>().createApplyMaskResult
	     (std::declval<unmasked_residual_type const &>() )
	     );

  using masked_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
      (std::declval<typename TrialSpaceType::basis_type const &>())
      );

public:
  // required aliases
  using state_type    = LspgStateType;
  using residual_type = masked_residual_type;
  using jacobian_type = masked_jac_action_result_type;

  LspgSteadyMaskedSystem() = delete;

  LspgSteadyMaskedSystem(TrialSpaceType & space,
			 const FomSystemType & fomSystem,
			 const MaskerType & masker)
    : space_(space),
      fomSystem_(fomSystem),
      fomState_(fomSystem.createState()),
      masker_(masker),
      unMaskedFomResidual_(fomSystem.createResidual()),
      unMaskedFomJacAction_(fomSystem.createApplyJacobianResult(space_.get().viewBasis()))
  {}

public:
  residual_type createResidual() const{
    auto tmp = fomSystem_.get().createResidual();
    return masker_.get().createApplyMaskResult(tmp);
  }

  jacobian_type createJacobian() const{
    auto tmp = fomSystem_.get().createApplyJacobianResult(space_.get().viewBasis());
    return masker_.get().createApplyMaskResult(tmp);
  }

  void residualAndJacobian(const state_type & lspgState,
			   residual_type & R,
			   jacobian_type & J,
			   bool recomputeJacobian = true) const
  {
    space_.get().mapFromReducedState(lspgState, fomState_);

    fomSystem_.get().residual(fomState_, unMaskedFomResidual_);
    masker_(unMaskedFomResidual_, R);

    if (recomputeJacobian){
      const auto & phi = space_.get().viewBasis();
      fomSystem_.get().applyJacobian(fomState_, phi, unMaskedFomJacAction_);
      masker_(unMaskedFomJacAction_, J);
    }
  }

protected:
  std::reference_wrapper<TrialSpaceType> space_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;

  std::reference_wrapper<const MaskerType> masker_;

  // UNMASKED fom R,J instances
  mutable unmasked_residual_type unMaskedFomResidual_;
  mutable unmasked_jac_action_result_type unMaskedFomJacAction_;
};

}}} // end pressio::rom::impl
#endif
