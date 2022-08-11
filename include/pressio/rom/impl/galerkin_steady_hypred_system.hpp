
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
  class TrialSpaceType,
  class FomSystemType,
  class HyperReductionOperator
  >
class GalerkinSteadyHypRedSystem
{

  using basis_type = typename TrialSpaceType::basis_type;

  // deduce from the fom object the type of result of
  // applying the Jacobian to the basis
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<basis_type const &>()) );

public:
  // required aliases
  using state_type    = ReducedStateType;
  using residual_type = ReducedResidualType;
  using jacobian_type = ReducedJacobianType;

  GalerkinSteadyHypRedSystem() = delete;

  GalerkinSteadyHypRedSystem(const TrialSpaceType & space,
			     const FomSystemType & fomSystem,
			     const HyperReductionOperator & hrOp)
    : space_(space),
      fomSystem_(fomSystem),
      hrOp_(hrOp),
      fomState_(fomSystem.createState()),
      fomResidual_(fomSystem.createResidual()),
      fomJacAction_(fomSystem.createApplyJacobianResult(space_.get().viewBasis()))
  {}

public:
  residual_type createResidual() const{
    // here we assume that the action of hyOp does not
    // produce a reduced residual different than the number of basis.
    // to be precise, we should compute: R = MJOP f(phi x)
    const auto & phi = space_.get().viewBasis();
    return impl::CreateGalerkinRhs<residual_type>()(phi);
  }

  jacobian_type createJacobian() const{
    // here we assume that the reduced jacobian is square matrix
    // defined by num of modes in basis.
    // to be precise, we should compute: J = MJOP df/dx(phi x) phi
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

    fomSystem_.get().residual(fomState_,  fomResidual_);
    hrOp_(fomResidual_, R);

    if (recomputeJacobian){
      fomSystem_.get().applyJacobian(fomState_, phi, fomJacAction_);
      hrOp_(fomJacAction_, J);
    }
  }

protected:
  std::reference_wrapper<const TrialSpaceType> space_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  std::reference_wrapper<const HyperReductionOperator> hrOp_;
  mutable typename FomSystemType::state_type fomState_;
  mutable typename FomSystemType::residual_type fomResidual_;
  mutable fom_jac_action_result_type fomJacAction_;
};

}}}
#endif
