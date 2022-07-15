
#ifndef PRESSIO_ROM_IMPL_GALERKIN_STEADY_DEFAULT_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_STEADY_DEFAULT_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

//
// phi^T fom_r(phi x) = 0
// so:
//   R = phi^T f_rom(phi x)
//   J = phi^T df/dx(phi x) phi
//
template <
  class ReducedStateType,
  class ResidualType,
  class JacobianType,
  class TrialSpaceType,
  class FomSystemType
  >
class GalerkinSteadyDefaultSystem
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

  GalerkinSteadyDefaultSystem() = delete;

  GalerkinSteadyDefaultSystem(const TrialSpaceType & space,
			      const FomSystemType & fomSystem)
    : space_(space), fomSystem_(fomSystem),
      fomState_(fomSystem.createState()),
      fomResidual_(fomSystem.createResidual()),
      fomJacAction_(fomSystem.createApplyJacobianResult(space_.get().viewBasis()))
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

    using scalar_t = typename ::pressio::Traits<basis_type>::scalar_type;
    using cnst = ::pressio::utils::Constants<scalar_t>;

    fomSystem_.get().residual(fomState_, fomResidual_);
    ::pressio::ops::product(::pressio::transpose(),
			    cnst::one(), phi, fomResidual_,
			    cnst::zero(), R);

    if (recomputeJacobian){
      fomSystem_.get().applyJacobian(fomState_, phi, fomJacAction_);
      ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
			      cnst::one(), phi, fomJacAction_,
			      cnst::zero(), J);
    }
  }

protected:
  std::reference_wrapper<const TrialSpaceType> space_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  mutable typename FomSystemType::residual_type fomResidual_;
  mutable fom_jac_action_result_type fomJacAction_;
};


}}} // end pressio::rom::impl
#endif
