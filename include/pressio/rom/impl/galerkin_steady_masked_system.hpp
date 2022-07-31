
#ifndef PRESSIO_ROM_IMPL_GALERKIN_STEADY_MASKED_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_STEADY_MASKED_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  class ReducedStateType,
  class ResidualType,
  class JacobianType,
  class TrialSpaceType,
  class FomSystemType,
  class ResidualMaskerType,
  class JacobianActionMaskerType,
  class HypRedOpType
  >
class GalerkinSteadyMaskedSystem
{

  using basis_type = typename TrialSpaceType::basis_type;
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<basis_type const &>()) );

public:
  // required aliases
  using state_type    = ReducedStateType;
  using residual_type = ResidualType;
  using jacobian_type = JacobianType;

  GalerkinSteadyMaskedSystem() = delete;
  GalerkinSteadyMaskedSystem(const TrialSpaceType & space,
			     const FomSystemType & fomSystem,
			     const ResidualMaskerType & rMasker,
			     const JacobianActionMaskerType & jaMasker,
			     const HypRedOpType & hrOp)
    : space_(space),
      fomSystem_(fomSystem),
      fomState_(fomSystem.createState()),
      rMasker_(rMasker),
      jaMasker_(jaMasker),
      hrOp_(hrOp),
      unMaskedFomResidual_(fomSystem.createResidual()),
      unMaskedFomJacAction_(fomSystem.createApplyJacobianResult(space_.get().viewBasis())),
      maskedFomResidual_(rMasker.createApplyMaskResult(unMaskedFomResidual_)),
      maskedFomJacAction_(jaMasker.createApplyMaskResult(unMaskedFomJacAction_))
  {}

public:
  residual_type createResidual() const{
    const auto & phi = space_.get().viewBasis();
    return impl::CreateGalerkinRhs<residual_type>()(phi);
  }

  jacobian_type createJacobian() const{
    const auto & phi = space_.get().viewBasis();
    return impl::CreateGalerkinJacobian<jacobian_type>()(phi);
  }

  void residualAndJacobian(const state_type & reducedState,
			   residual_type & R,
			   jacobian_type & J,
			   bool recomputeJacobian = true) const
  {
    const auto & phi = space_.get().viewBasis();
    space_.get().mapFromReducedState(reducedState, fomState_);

    fomSystem_.get().residual(fomState_, unMaskedFomResidual_);
    rMasker_(unMaskedFomResidual_, maskedFomResidual_);
    hrOp_(maskedFomResidual_, R);

    if (recomputeJacobian){
      fomSystem_.get().applyJacobian(fomState_, phi, unMaskedFomJacAction_);
      jaMasker_(unMaskedFomJacAction_, maskedFomJacAction_);
      hrOp_(maskedFomJacAction_, J);
    }
  }

private:
  std::reference_wrapper<const TrialSpaceType> space_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;

  // masker and hyp-red operator
  std::reference_wrapper<const ResidualMaskerType> rMasker_;
  std::reference_wrapper<const JacobianActionMaskerType> jaMasker_;
  std::reference_wrapper<const HypRedOpType> hrOp_;

  // UNMASKED fom R,J instances
  mutable typename FomSystemType::residual_type unMaskedFomResidual_;
  mutable fom_jac_action_result_type unMaskedFomJacAction_;

  // MASKED fom R,J instances
  mutable typename FomSystemType::residual_type maskedFomResidual_;
  mutable fom_jac_action_result_type maskedFomJacAction_;
};

}}}
#endif
