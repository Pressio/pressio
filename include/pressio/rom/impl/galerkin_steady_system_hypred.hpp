
#ifndef PRESSIO_ROM_IMPL_GALERKIN_STEADY_HYPER_REDUCED_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_STEADY_HYPER_REDUCED_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
hypred Galerkin problem represents:

   hypredOp fom_r(phi x) = 0

- fom_r is the fom "residual"
- phi is the basis

From this we get a "reduced" residual/jacobian:
R = phi^T hypredOp fom_r(phi x)
J = phi^T hypredOp dfom_r/dx(phi x) phi
*/
template <
  class ReducedStateType,
  class ReducedResidualType,
  class ReducedJacobianType,
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReductionOperator
  >
class GalerkinSteadyHypRedSystem
{

  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;

  // deduce from the fom object the type of result of
  // applying the Jacobian to the basis
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<basis_matrix_type const &>()) );

public:
  // required aliases
  using state_type    = ReducedStateType;
  using residual_type = ReducedResidualType;
  using jacobian_type = ReducedJacobianType;

  GalerkinSteadyHypRedSystem() = delete;

  GalerkinSteadyHypRedSystem(const TrialSubspaceType & trialSubspace,
			     const FomSystemType & fomSystem,
			     const HyperReductionOperator & hrOp)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      hrOp_(hrOp),
      fomResidual_(fomSystem.createResidual()),
      fomJacAction_(fomSystem.createApplyJacobianResult(trialSubspace_.get().basisOfTranslatedSpace()))
  {}

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  residual_type createResidual() const{
    // here we assume that the action of hyOp does not
    // produce a reduced residual different than the number of basis.
    // to be precise, we should compute: R = MJOP f(phi x)
    return impl::CreateGalerkinRhs<residual_type>()(trialSubspace_.get().dimension());
  }

  jacobian_type createJacobian() const{
    // here we assume that the reduced jacobian is square matrix
    // defined by num of modes in basis.
    // to be precise, we should compute: J = MJOP df/dx(phi x) phi
    return impl::CreateGalerkinJacobian<jacobian_type>()(trialSubspace_.get().dimension());
  }

  void residualAndJacobian(const state_type & reducedState,
			   residual_type & reducedResidual,
			   jacobian_type & reducedJacobian,
			   bool computeJacobian) const
  {

    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    fomSystem_.get().residual(fomState_,  fomResidual_);
    hrOp_(fomResidual_, reducedResidual);

    if (computeJacobian){
      fomSystem_.get().applyJacobian(fomState_, phi, fomJacAction_);
      hrOp_(fomJacAction_, reducedJacobian);
    }
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReductionOperator> hrOp_;
  mutable typename FomSystemType::residual_type fomResidual_;
  mutable fom_jac_action_result_type fomJacAction_;
};

}}}
#endif
