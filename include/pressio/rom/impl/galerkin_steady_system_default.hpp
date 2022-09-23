
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
  class TrialSubspaceType,
  class FomSystemType
  >
class GalerkinSteadyDefaultSystem
{

  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;

  // deduce from the fom object the type of result of
  // applying the Jacobian to the basis
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<basis_matrix_type const &>()) );

public:
  // aliases required by the pressio solvers
  using state_type    = ReducedStateType;
  using residual_type = ReducedResidualType;
  using jacobian_type = ReducedJacobianType;

  GalerkinSteadyDefaultSystem() = delete;

  GalerkinSteadyDefaultSystem(const TrialSubspaceType & trialSubspace,
			      const FomSystemType & fomSystem)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      fomResidual_(fomSystem.createResidual()),
      fomJacAction_(fomSystem.createApplyJacobianResult(trialSubspace_.get().basisOfTranslatedSpace()))
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

    using phi_scalar_t = typename ::pressio::Traits<basis_matrix_type>::scalar_type;
    constexpr auto alpha = ::pressio::utils::Constants<phi_scalar_t>::one();
    using R_scalar_t = typename ::pressio::Traits<residual_type>::scalar_type;
    constexpr auto beta = ::pressio::utils::Constants<R_scalar_t>::zero();
    fomSystem_.get().residual(fomState_, fomResidual_);
    ::pressio::ops::product(::pressio::transpose(),
			    alpha, phi, fomResidual_, beta, reducedResidual);

    if (computeJacobian){
      fomSystem_.get().applyJacobian(fomState_, phi, fomJacAction_);

      using J_scalar_t = typename ::pressio::Traits<jacobian_type>::scalar_type;
      constexpr auto beta = ::pressio::utils::Constants<J_scalar_t>::zero();
      ::pressio::ops::product(::pressio::transpose(),
			      ::pressio::nontranspose(),
			      alpha, phi, fomJacAction_, beta, reducedJacobian);
    }
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  mutable typename FomSystemType::residual_type fomResidual_;
  mutable fom_jac_action_result_type fomJacAction_;
};

}}} // end pressio::rom::impl
#endif
