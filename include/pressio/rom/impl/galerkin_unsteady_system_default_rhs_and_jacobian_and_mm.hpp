
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_DEFAULT_RHS_AND_JACOBIAN_AND_MM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_DEFAULT_RHS_AND_JACOBIAN_AND_MM_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  default implicit galerkin system represents:

     phi^T M(y,t,...) phi d hat{y}/dt = phi^T fom_rhs(phi*hat{y}, ...)

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class ReducedJacobianType,
  class ReducedMassMatrixType,
  class TrialSubspaceType,
  class FomSystemType
  >
class GalerkinDefaultOdeSystemRhsJacobianMassMatrix
{

  // deduce from the fom object the type of result of
  // applying the Jacobian to the basis
  using basis_matrix_type = typename TrialSubspaceType::basis_matrix_type;
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createResultOfJacobianActionOn
	     (std::declval<basis_matrix_type const &>()) );

  using fom_mm_action_result_type =
    decltype(std::declval<FomSystemType const>().createResultOfMassMatrixActionOn
	     (std::declval<basis_matrix_type const &>()) );

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using rhs_type      = ReducedRhsType;
  using jacobian_type             = ReducedJacobianType;
  using mass_matrix_type	  = ReducedMassMatrixType;

  GalerkinDefaultOdeSystemRhsJacobianMassMatrix(const TrialSubspaceType & trialSubspace,
						const FomSystemType & fomSystem)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState()),
      fomRhs_(fomSystem.createRhs()),
      fomJacAction_(fomSystem.createResultOfJacobianActionOn(trialSubspace_.get().basisOfTranslatedSpace())),
      fomMMAction_(fomSystem.createResultOfMassMatrixActionOn(trialSubspace_.get().basisOfTranslatedSpace()))
  {}

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  rhs_type createRhs() const{
    return impl::CreateGalerkinRhs<rhs_type>()(trialSubspace_.get().dimension());
  }

  jacobian_type createJacobian() const{
    return impl::CreateGalerkinJacobian<jacobian_type>()(trialSubspace_.get().dimension());
  }

  mass_matrix_type createMassMatrix() const{
    return impl::CreateGalerkinMassMatrix<mass_matrix_type>()(trialSubspace_.get().dimension());
  }

  void massMatrixAndRhsAndJacobian(const state_type & reducedState,
				   const IndVarType & rhsEvaluationTime,
				   mass_matrix_type & reducedMassMatrix,
				   rhs_type & reducedRhs,
				   std::optional<jacobian_type*> reducedJacobian) const
  {

    // reconstruct fom state fomState = phi*reducedState
    trialSubspace_.get().mapFromReducedState(reducedState, fomState_);

    // evaluate fomRhs
    fomSystem_.get().rhs(fomState_, rhsEvaluationTime, fomRhs_);

    // compute the reduced rhs
    const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
    using phi_scalar_t = typename ::pressio::Traits<basis_matrix_type>::scalar_type;
    constexpr auto alpha = static_cast<phi_scalar_t>(1);
    using rhs_scalar_t = typename ::pressio::Traits<rhs_type>::scalar_type;
    constexpr auto beta = static_cast<rhs_scalar_t>(0);
    ::pressio::ops::product(::pressio::transpose(),
			    alpha, phi, fomRhs_,
			    beta, reducedRhs);

    // compute the reduced mass matrix
    fomSystem_.get().applyMassMatrix(fomState_, phi, rhsEvaluationTime, fomMMAction_);
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
			    alpha, phi, fomMMAction_,
			    beta, reducedMassMatrix);

    if (reducedJacobian){
      // evaluate fom jacobian action: fomJacAction_ = fom_J * phi
      fomSystem_.get().applyJacobian(fomState_, phi, rhsEvaluationTime, fomJacAction_);

      // compute the reduced jacobian
      constexpr auto alpha = static_cast<phi_scalar_t>(1);
      constexpr auto beta = static_cast<rhs_scalar_t>(0);
      ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
			      alpha, phi, fomJacAction_,
			      beta,
			      *reducedJacobian.value());
    }
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  mutable typename FomSystemType::rhs_type fomRhs_;
  mutable fom_jac_action_result_type fomJacAction_;
  mutable fom_mm_action_result_type fomMMAction_;
};

}}} // end pressio::rom::impl
#endif  // PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEM_DEFAULT_RHS_AND_JACOBIAN_AND_MM_HPP_
