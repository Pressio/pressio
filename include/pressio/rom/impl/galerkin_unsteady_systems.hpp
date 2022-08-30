
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEMS_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEMS_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
  default explicit galerkin system represents:

     d hat{y}/dt = phi^T fom_rhs(phi*hat{y}, ...)

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis

*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class TrialSpaceType,
  class FomSystemType
  >
class GalerkinDefaultOdeSystemOnlyRhs
{
  using basis_type = typename TrialSpaceType::basis_type;

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using right_hand_side_type      = ReducedRhsType;

  GalerkinDefaultOdeSystemOnlyRhs() = delete;

  GalerkinDefaultOdeSystemOnlyRhs(const TrialSpaceType & trialSpace,
				  const FomSystemType & fomSystem)
    : trialSpace_(trialSpace),
      fomSystem_(fomSystem),
      fomState_(fomSystem.createState()),
      fomRhs_(fomSystem.createRightHandSide())
  {}

public:
  state_type createState() const{
    return trialSpace_.get().createReducedState();
  }

  right_hand_side_type createRightHandSide() const{
    const auto & phi = trialSpace_.get().viewBasis();
    return impl::CreateGalerkinRhs<right_hand_side_type>()(phi);
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

protected:
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  mutable typename FomSystemType::right_hand_side_type fomRhs_;
};


/*
  hyp-red explicit galerkin system represents:

     d hat{y}/dt = HrOp fom_rhs(phi*hat{y}, ...)

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
- HrOp is the hyper-red operator
*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class TrialSpaceType,
  class FomSystemType,
  class HyperReductionOperator
  >
class GalerkinHyperReducedOdeSystemOnlyRhs
{

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using right_hand_side_type      = ReducedRhsType;

  GalerkinHyperReducedOdeSystemOnlyRhs() = delete;

  GalerkinHyperReducedOdeSystemOnlyRhs(const TrialSpaceType & trialSpace,
				       const FomSystemType & fomSystem,
				       const HyperReductionOperator & hrOp)
    : trialSpace_(trialSpace),
      fomSystem_(fomSystem),
      hrOp_(hrOp),
      fomState_(fomSystem.createState()),
      fomRhs_(fomSystem.createRightHandSide())
  {}

public:
  state_type createState() const{
    return trialSpace_.get().createReducedState();
  }

  right_hand_side_type createRightHandSide() const{
    const auto & phi = trialSpace_.get().viewBasis();
    return impl::CreateGalerkinRhs<right_hand_side_type>()(phi);
  }

  void rightHandSide(const state_type & reducedState,
		     const IndVarType & rhsEvaluationTime,
		     right_hand_side_type & reducedRhs) const
  {

    // reconstruct fom state fomState = phi*reducedState
    trialSpace_.get().mapFromReducedState(reducedState, fomState_);

    // evaluate fomRhs
    fomSystem_.get().rightHandSide(fomState_, rhsEvaluationTime, fomRhs_);

    // evaluate reduced rhs
    hrOp_(fomRhs_, rhsEvaluationTime, reducedRhs);
  }

protected:
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  std::reference_wrapper<const HyperReductionOperator> hrOp_;
  mutable typename FomSystemType::state_type fomState_;
  mutable typename FomSystemType::right_hand_side_type fomRhs_;
};


/*
  masked explicit galerkin system represents:

     d hat{y}/dt = HrOp mask( fom_rhs(phi*hat{y}, ...) )

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
- HrOp is the hyper-red operator
*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class TrialSpaceType,
  class FomSystemType,
  class RhsMaskerType,
  class HyperReductionOperator
  >
class GalerkinMaskedOdeSystemOnlyRhs
{
  // deduce types
  using unmasked_fom_rhs_type = typename FomSystemType::right_hand_side_type;
  using masked_fom_rhs_type = typename RhsMaskerType::result_type;

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using right_hand_side_type      = ReducedRhsType;

  GalerkinMaskedOdeSystemOnlyRhs() = delete;

  GalerkinMaskedOdeSystemOnlyRhs(const TrialSpaceType & trialSpace,
				 const FomSystemType & fomSystem,
				 const RhsMaskerType & rhsMasker,
				 const HyperReductionOperator & hrOp)
    : trialSpace_(trialSpace),
      fomSystem_(fomSystem),
      fomState_(fomSystem.createState()),
      hrOp_(hrOp),
      rhsMasker_(rhsMasker),
      unMaskedFomRhs_(fomSystem.createRightHandSide()),
      maskedFomRhs_(rhsMasker.createApplyMaskResult(unMaskedFomRhs_))
  {}

public:
  state_type createState() const{
    return trialSpace_.get().createReducedState();
  }

  right_hand_side_type createRightHandSide() const{
    const auto & phi = trialSpace_.get().viewBasis();
    return impl::CreateGalerkinRhs<right_hand_side_type>()(phi);
  }

  void rightHandSide(const state_type & reducedState,
		     const IndVarType & rhsEvaluationTime,
		     right_hand_side_type & reducedRhs) const
  {

    // reconstruct fom state fomState = phi*reducedState
    trialSpace_.get().mapFromReducedState(reducedState, fomState_);

    // evaluate fomRhs and mask it
    fomSystem_.get().rightHandSide(fomState_, rhsEvaluationTime, unMaskedFomRhs_);
    rhsMasker_(unMaskedFomRhs_, maskedFomRhs_);

    // evaluate reduced rhs
    hrOp_(maskedFomRhs_, rhsEvaluationTime, reducedRhs);
  }

protected:
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReductionOperator> hrOp_;

  std::reference_wrapper<const RhsMaskerType> rhsMasker_;
  mutable unmasked_fom_rhs_type unMaskedFomRhs_;
  mutable masked_fom_rhs_type maskedFomRhs_;
};


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
      fomState_(fomSystem.createState()),
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



/*
  hypred implicit galerkin system represents:

     d hat{y}/dt = hrOp fom_rhs(phi*hat{y}, ...)

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
- hrOp is the hypred operator
so that it boils down to:

rhs = hrOp fom_rhs(phi*hat{y}, ...)
rhs_jacobian = hrOp d(fom_rhs(phi*hat{y}, ...))/dy phi

*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class ReducedJacobianType,
  class TrialSpaceType,
  class FomSystemType,
  class HyperReductionOperator
  >
class GalerkinHypRedOdeSystemRhsAndJacobian
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

  GalerkinHypRedOdeSystemRhsAndJacobian() = delete;

  GalerkinHypRedOdeSystemRhsAndJacobian(const TrialSpaceType & trialSpace,
					const FomSystemType & fomSystem,
					const HyperReductionOperator & hrOp)
    : trialSpace_(trialSpace),
      fomSystem_(fomSystem),
      hrOp_(hrOp),
      fomState_(fomSystem.createState()),
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

    // evaluate fomRhs and apply hrOp
    fomSystem_.get().rightHandSide(fomState_, rhsEvaluationTime, fomRhs_);

    // compute the reduced rhs
    hrOp_(fomRhs_, rhsEvaluationTime, reducedRhs);
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
    hrOp_(fomJacAction_, rhsEvaluationTime, reducedJacobian);
  }

protected:
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  std::reference_wrapper<const HyperReductionOperator> hrOp_;
  mutable typename FomSystemType::state_type fomState_;
  mutable typename FomSystemType::right_hand_side_type fomRhs_;
  mutable fom_jac_action_result_type fomJacAction_;
};


/*
  masked implicit galerkin system represents:

     d hat{y}/dt = hrOp masker ( fom_rhs(phi*hat{y}, ...) )

- hat{y} is the reduced state
- fom_rhs is the fom RHS
- phi is the basis
- hrOp is the hypred operator
so that it boils down to:

rhs = hrOp masked(fom_rhs(phi*hat{y}, ...))
rhs_jacobian = hrOp masked(d(fom_rhs(phi*hat{y}, ...))/dy phi)

*/
template <
  class IndVarType,
  class ReducedStateType,
  class ReducedRhsType,
  class ReducedJacobianType,
  class TrialSpaceType,
  class FomSystemType,
  class RhsMaskerType,
  class JacobianActionMaskerType,
  class HyperReductionOperator
  >
class GalerkinMaskedOdeSystemRhsAndJacobian
{
  using basis_type = typename TrialSpaceType::basis_type;

  // deduce the unmasked types
  using unmasked_fom_rhs_type = typename FomSystemType::right_hand_side_type;
  using unmasked_fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<basis_type const &>()));

  // deduce the masked types
  using masked_fom_rhs_type = typename RhsMaskerType::result_type;
  using masked_fom_jac_action_result_type = typename JacobianActionMaskerType::result_type;

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using right_hand_side_type      = ReducedRhsType;
  using jacobian_type = ReducedJacobianType;

  GalerkinMaskedOdeSystemRhsAndJacobian() = delete;

  GalerkinMaskedOdeSystemRhsAndJacobian(const TrialSpaceType & trialSpace,
					const FomSystemType & fomSystem,
					const RhsMaskerType & rhsMasker,
					const JacobianActionMaskerType & jaMasker,
					const HyperReductionOperator & hrOp)
    : trialSpace_(trialSpace),
      fomSystem_(fomSystem),
      fomState_(fomSystem.createState()),
      hrOp_(hrOp),
      rhsMasker_(rhsMasker),
      jaMasker_(jaMasker),
      unMaskedFomRhs_(fomSystem.createRightHandSide()),
      unMaskedFomJacAction_(fomSystem.createApplyJacobianResult(trialSpace_.get().viewBasis())),
      maskedFomRhs_(rhsMasker.createApplyMaskResult(unMaskedFomRhs_)),
      maskedFomJacAction_(jaMasker.createApplyMaskResult(unMaskedFomJacAction_))
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

    trialSpace_.get().mapFromReducedState(reducedState, fomState_);
    fomSystem_.get().rightHandSide(fomState_, rhsEvaluationTime, unMaskedFomRhs_);
    rhsMasker_(unMaskedFomRhs_, maskedFomRhs_);
    hrOp_(maskedFomRhs_, rhsEvaluationTime, reducedRhs);
  }

  void jacobian(const state_type & reducedState,
		const IndVarType & rhsEvaluationTime,
		jacobian_type & reducedJacobian) const

  {
    trialSpace_.get().mapFromReducedState(reducedState, fomState_);

    const auto & phi = trialSpace_.get().viewBasis();
    fomSystem_.get().applyJacobian(fomState_, phi, rhsEvaluationTime, unMaskedFomJacAction_);
    jaMasker_(unMaskedFomJacAction_, maskedFomJacAction_);
    hrOp_(maskedFomJacAction_, rhsEvaluationTime, reducedJacobian);
  }

protected:
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  std::reference_wrapper<const HyperReductionOperator> hrOp_;

  // masker
  std::reference_wrapper<const RhsMaskerType> rhsMasker_;
  std::reference_wrapper<const JacobianActionMaskerType> jaMasker_;

  // UNMASKED objects
  mutable unmasked_fom_rhs_type unMaskedFomRhs_;
  mutable unmasked_fom_jac_action_result_type unMaskedFomJacAction_;

  // MASKED objects
  mutable masked_fom_rhs_type maskedFomRhs_;
  mutable masked_fom_jac_action_result_type maskedFomJacAction_;
};

}}} // end pressio::rom::impl
#endif
