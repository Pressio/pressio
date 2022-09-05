
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_DEFAULT_RHS_AND_JACOBIAN_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_DEFAULT_RHS_AND_JACOBIAN_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  default implicit galerkin system represents:

     d hat{y}/dt = phi^T fom_rhs(phi*hat{y}, ...)

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
so that it boils down to:

rhs = phi^T fom_rhs(phi*hat{y}, ...)
rhs_jacobian = phi^T d(fom_rhs(phi*hat{y}, ...))/dy phi

*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class ReducedJacobianType,
  class TrialSpaceType,
  class FomSystemType
  >
class GalerkinDefaultOdeSystemRhsAndJacobian
{

  // deduce from the fom object the type of result of
  // applying the Jacobian to the basis
  using basis_type = typename TrialSpaceType::basis_type;
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<basis_type const &>()) );

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using right_hand_side_type      = ReducedRhsType;
  using jacobian_type = ReducedJacobianType;

  GalerkinDefaultOdeSystemRhsAndJacobian() = delete;

  GalerkinDefaultOdeSystemRhsAndJacobian(const TrialSpaceType & trialSpace,
					 const FomSystemType & fomSystem)
    : trialSpace_(trialSpace),
      fomSystem_(fomSystem),
      fomState_(trialSpace.createFullState()),
      fomRhs_(fomSystem.createRightHandSide()),
      fomJacAction_(fomSystem.createApplyJacobianResult(trialSpace_.get().viewBasis()))
  {}

public:
  state_type createState() const{
    return trialSpace_.get().createReducedState();
  }

  right_hand_side_type createRightHandSide() const{
    const auto & phi = trialSpace_.get().viewBasis();
    return impl::CreateGalerkinRhs<right_hand_side_type>()(phi);
  }

  jacobian_type createJacobian() const{
    const auto & phi = trialSpace_.get().viewBasis();
    return impl::CreateGalerkinJacobian<jacobian_type>()(phi);
  }

  void rightHandSide(const state_type & reducedState,
		     const IndVarType & rhsEvaluationTime,
		     right_hand_side_type & reducedRhs) const
  {

    // reconstruct fom state fomState = phi*reducedState
    trialSpace_.get().mapFromReducedState(reducedState, fomState_);

    // evaluate fomRhs
    fomSystem_.get().rightHandSide(fomState_, rhsEvaluationTime, fomRhs_);

    // compute the reduced rhs
    const auto & phi = trialSpace_.get().viewBasis();
    using phi_scalar_t = typename ::pressio::Traits<basis_type>::scalar_type;
    constexpr auto alpha = ::pressio::utils::Constants<phi_scalar_t>::one();
    using rhs_scalar_t = typename ::pressio::Traits<right_hand_side_type>::scalar_type;
    constexpr auto beta = ::pressio::utils::Constants<rhs_scalar_t>::zero();
    ::pressio::ops::product(::pressio::transpose(),
			    alpha, phi, fomRhs_,
			    beta, reducedRhs);
  }

  void jacobian(const state_type & reducedState,
		const IndVarType & rhsEvaluationTime,
		jacobian_type & reducedJacobian) const

  {
    // reconstruct fom state fomState = phi*reducedState
    trialSpace_.get().mapFromReducedState(reducedState, fomState_);

    const auto & phi = trialSpace_.get().viewBasis();

    // evaluate fom jacobian action: fomJacAction_ = fom_J * phi
    fomSystem_.get().applyJacobian(fomState_, phi, rhsEvaluationTime, fomJacAction_);

    // compute the reduced jacobian
    using phi_scalar_t = typename ::pressio::Traits<basis_type>::scalar_type;
    constexpr auto alpha = ::pressio::utils::Constants<phi_scalar_t>::one();
    using rhs_scalar_t = typename ::pressio::Traits<right_hand_side_type>::scalar_type;
    constexpr auto beta = ::pressio::utils::Constants<rhs_scalar_t>::zero();
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
			    alpha, phi, fomJacAction_,
			    beta, reducedJacobian);
  }

protected:
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  mutable typename FomSystemType::right_hand_side_type fomRhs_;
  mutable fom_jac_action_result_type fomJacAction_;
};

}}} // end pressio::rom::impl
#endif
