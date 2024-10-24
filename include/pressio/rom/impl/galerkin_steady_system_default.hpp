
#ifndef PRESSIO_ROM_IMPL_GALERKIN_STEADY_SYSTEM_DEFAULT_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_STEADY_SYSTEM_DEFAULT_HPP_

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
    decltype(std::declval<FomSystemType const>().createResultOfJacobianActionOn
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
      fomJacAction_(fomSystem.createResultOfJacobianActionOn(trialSubspace_.get().basisOfTranslatedSpace()))
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
#ifdef PRESSIO_ENABLE_CXX17
			   std::optional<jacobian_type*> reducedJacobian) const
#else
                           jacobian_type* reducedJacobian) const
#endif
  {

    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    using phi_scalar_t = typename ::pressio::Traits<basis_matrix_type>::scalar_type;
    constexpr auto alpha = static_cast<phi_scalar_t>(1);
    using R_scalar_t = typename ::pressio::Traits<residual_type>::scalar_type;
    constexpr auto beta = static_cast<R_scalar_t>(0);

#ifdef PRESSIO_ENABLE_CXX17
    fomSystem_.get().residualAndJacobianAction(fomState_, fomResidual_, phi,
					       std::optional<fom_jac_action_result_type *>(&fomJacAction_));
#else
    fomSystem_.get().residualAndJacobianAction(fomState_, fomResidual_, phi, &fomJacAction_);
#endif

    ::pressio::ops::product(::pressio::transpose(),
			    alpha, phi, fomResidual_, beta, reducedResidual);

    if (reducedJacobian){
      using J_scalar_t = typename ::pressio::Traits<jacobian_type>::scalar_type;
      constexpr auto beta = static_cast<J_scalar_t>(0);
      ::pressio::ops::product(::pressio::transpose(),
			      ::pressio::nontranspose(),
			      alpha, phi, fomJacAction_, beta,
#ifdef PRESSIO_ENABLE_CXX17
			      *reducedJacobian.value());
#else
                              *reducedJacobian);
#endif
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
#endif  // PRESSIO_ROM_IMPL_GALERKIN_STEADY_SYSTEM_DEFAULT_HPP_
