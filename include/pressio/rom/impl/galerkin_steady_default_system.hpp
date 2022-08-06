
#ifndef PRESSIO_ROM_IMPL_GALERKIN_STEADY_DEFAULT_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_STEADY_DEFAULT_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
default Galerkin problem represents:

  phi^T fom_r(phi x) = 0

- fom_r is the fom "residual"
- phi is the basis

From this we get a "reduced" residual/jacobian:
R = phi^T fom_r(phi x)
J = phi^T dfom_r/dx(phi x) phi
*/
template <
  class ReducedStateType,
  class ReducedResidualType,
  class ReducedJacobianType,
  class TrialSpaceType,
  class FomSystemType
  >
class GalerkinSteadyDefaultSystem
{

  using basis_type = typename TrialSpaceType::basis_type;

  // deduce from the fom object the type of result of
  // applying the Jacobian to the basis
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<basis_type const &>()) );

public:
  // aliases required by the pressio solvers
  using state_type    = ReducedStateType;
  using residual_type = ReducedResidualType;
  using jacobian_type = ReducedJacobianType;

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

    using phi_scalar_t = typename ::pressio::Traits<basis_type>::scalar_type;
    constexpr auto alpha = ::pressio::utils::Constants<phi_scalar_t>::one();
    using R_scalar_t = typename ::pressio::Traits<residual_type>::scalar_type;
    constexpr auto beta = ::pressio::utils::Constants<R_scalar_t>::zero();
    fomSystem_.get().residual(fomState_, fomResidual_);
    ::pressio::ops::product(::pressio::transpose(),
			    alpha, phi, fomResidual_, beta, R);

    if (recomputeJacobian){
      fomSystem_.get().applyJacobian(fomState_, phi, fomJacAction_);

      using J_scalar_t = typename ::pressio::Traits<jacobian_type>::scalar_type;
      constexpr auto beta = ::pressio::utils::Constants<J_scalar_t>::zero();
      ::pressio::ops::product(::pressio::transpose(),
			      ::pressio::nontranspose(),
			      alpha, phi, fomJacAction_, beta, J);
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
