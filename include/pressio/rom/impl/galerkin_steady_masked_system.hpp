
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
  class TrialSpaceType,
  class FomSystemType,
  class ResidualMaskerType,
  class JacobianActionMaskerType,
  class HypRedOpType
  >
class GalerkinSteadyMaskedSystem
{

  using basis_type = typename TrialSpaceType::basis_type;

  // deduce the unmasked types
  using unmasked_fom_residual_type = typename FomSystemType::residual_type;
  using unmasked_fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<basis_type const &>()));

  // deduce the masked types
  using masked_fom_residual_type = typename ResidualMaskerType::result_type;
  using masked_fom_jac_action_result_type = typename JacobianActionMaskerType::result_type;

public:
  // required aliases
  using state_type    = ReducedStateType;
  using residual_type = ReducedResidualType;
  using jacobian_type = ReducedJacobianType;

  GalerkinSteadyMaskedSystem() = delete;
  GalerkinSteadyMaskedSystem(const TrialSpaceType & trialSpace,
			     const FomSystemType & fomSystem,
			     const ResidualMaskerType & rMasker,
			     const JacobianActionMaskerType & jaMasker,
			     const HypRedOpType & hrOp)
    : trialSpace_(trialSpace),
      fomSystem_(fomSystem),
      fomState_(trialSpace.createFullState()),
      hrOp_(hrOp),
      rMasker_(rMasker),
      jaMasker_(jaMasker),
      unMaskedFomResidual_(fomSystem.createResidual()),
      unMaskedFomJacAction_(fomSystem.createApplyJacobianResult(trialSpace_.get().viewBasis())),
      maskedFomResidual_(rMasker.createApplyMaskResult(unMaskedFomResidual_)),
      maskedFomJacAction_(jaMasker.createApplyMaskResult(unMaskedFomJacAction_))
  {}

public:
  state_type createState() const{
    return trialSpace_.get().createReducedState();
  }

  residual_type createResidual() const{
    const auto & phi = trialSpace_.get().viewBasis();
    return impl::CreateGalerkinRhs<residual_type>()(phi);
  }

  jacobian_type createJacobian() const{
    const auto & phi = trialSpace_.get().viewBasis();
    return impl::CreateGalerkinJacobian<jacobian_type>()(phi);
  }

  void residualAndJacobian(const state_type & reducedState,
			   residual_type & reducedResidual,
			   jacobian_type & reducedJacobian,
			   bool computeJacobian) const
  {
    const auto & phi = trialSpace_.get().viewBasis();
    trialSpace_.get().mapFromReducedState(reducedState, fomState_);

    fomSystem_.get().residual(fomState_, unMaskedFomResidual_);
    rMasker_(unMaskedFomResidual_, maskedFomResidual_);
    hrOp_(maskedFomResidual_, reducedResidual);

    if (computeJacobian){
      fomSystem_.get().applyJacobian(fomState_, phi, unMaskedFomJacAction_);
      jaMasker_(unMaskedFomJacAction_, maskedFomJacAction_);
      hrOp_(maskedFomJacAction_, reducedJacobian);
    }
  }

private:
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HypRedOpType> hrOp_;

  // masker
  std::reference_wrapper<const ResidualMaskerType> rMasker_;
  std::reference_wrapper<const JacobianActionMaskerType> jaMasker_;

  // UNMASKED fom R,J instances
  mutable unmasked_fom_residual_type unMaskedFomResidual_;
  mutable unmasked_fom_jac_action_result_type unMaskedFomJacAction_;

  // MASKED fom R,J instances
  mutable masked_fom_residual_type maskedFomResidual_;
  mutable masked_fom_jac_action_result_type maskedFomJacAction_;
};

}}}
#endif
