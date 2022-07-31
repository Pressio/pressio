
#ifndef PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEMS_HPP_
#define PRESSIO_ROM_IMPL_GALERKIN_UNSTEADY_SYSTEMS_HPP_

namespace pressio{ namespace rom{ namespace impl{

//
// system implementing:
// d hat{y}/dt = phi^T f
//
template <
  class IndVarType,
  class ReducedStateType,
  class GalerkinRhsType,
  class TrialSpaceType,
  class FomSystemType
  >
class GalerkinUnsteadyDefaultOdeSystem
{
public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using right_hand_side_type      = GalerkinRhsType;

  GalerkinUnsteadyDefaultOdeSystem() = delete;

  GalerkinUnsteadyDefaultOdeSystem(const TrialSpaceType & space,
				   const FomSystemType & fomSystem)
    : space_(space),
      fomSystem_(fomSystem),
      fomState_(fomSystem.createState()),
      fomRhs_(fomSystem.createRightHandSide())
  {}

public:
  state_type createState() const{
    return space_.get().createReducedState();
  }

  right_hand_side_type createRightHandSide() const{
    const auto & phi = space_.get().viewBasis();
    return impl::CreateGalerkinRhs<right_hand_side_type>()(phi);
  }

  void rightHandSide(const state_type & galerkinState,
		     const IndVarType & timeForRhsEvaluation,
		     right_hand_side_type & galerkinRhs) const
  {
    space_.get().mapFromReducedState(galerkinState, fomState_);
    fomSystem_.get().rightHandSide(fomState_, timeForRhsEvaluation, fomRhs_);

    const auto & phi = space_.get().viewBasis();

    using phi_type = typename TrialSpaceType::manifold_jacobian_type;
    using scalar_t = typename ::pressio::Traits<phi_type>::scalar_type;
    using cnst = ::pressio::utils::Constants<scalar_t>;
    ::pressio::ops::product(::pressio::transpose(), cnst::one(),
			    phi, fomRhs_, cnst::zero(), galerkinRhs);
  }

protected:
  std::reference_wrapper<const TrialSpaceType> space_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
  mutable typename FomSystemType::right_hand_side_type fomRhs_;
};


//
// system class implementing:
// phi^t M(y) phi d hat{y}/dt = phi^T f
//
template <
  class IndVarType,
  class ReducedStateType,
  class GalerkinRhsType,
  class GalerkinMassMatrixType,
  class TrialSpaceType,
  class FomSystemType
  >
class GalerkinUnsteadyDefaultOdeSystemWithMassMatrix
  : public GalerkinUnsteadyDefaultOdeSystem<
  IndVarType, ReducedStateType, GalerkinRhsType, TrialSpaceType, FomSystemType>
{
  using base_t = GalerkinUnsteadyDefaultOdeSystem<
    IndVarType, ReducedStateType, GalerkinRhsType,
    TrialSpaceType, FomSystemType>;

  // deduce the return type of the action of the mass matrix
  using basis_type = typename TrialSpaceType::basis_type;
  using fom_mm_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyMassMatrixResult
	     (std::declval<basis_type const &>())
	     );

public:
  // required aliases
  using typename base_t::independent_variable_type;
  using typename base_t::state_type;
  using typename base_t::right_hand_side_type;
  using mass_matrix_type  = GalerkinMassMatrixType;

  GalerkinUnsteadyDefaultOdeSystemWithMassMatrix() = delete;

  GalerkinUnsteadyDefaultOdeSystemWithMassMatrix(const TrialSpaceType & space,
						 const FomSystemType & fomSystem)
    : base_t(space, fomSystem),
      fomMMAction_(fomSystem.createApplyMassMatrixResult
		   (space_.get().viewBasis())
		   )
  {}

public:
  mass_matrix_type createMassMatrix() const{
    const auto & phi = space_.get().viewBasis();
    return impl::CreateGalerkinMassMatrix<mass_matrix_type>()(phi);
  }

  void massMatrix(const state_type & galerkinState,
		  const independent_variable_type & timeForEvaluation,
		  mass_matrix_type & galerkinMM) const
  {

    space_.get().mapFromReducedState(galerkinState, fomState_);

    const auto & phi = space_.get().viewBasis();
    fomSystem_.get().applyMassMatrix(fomState_, phi, timeForEvaluation, fomMMAction_);

    using phi_type = typename TrialSpaceType::basis_type;
    using scalar_t = typename ::pressio::Traits<phi_type>::scalar_type;
    using cnst = ::pressio::utils::Constants<scalar_t>;
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
			    cnst::one(), phi, fomMMAction_, cnst::zero(), galerkinMM);
  }

private:
  using base_t::space_;
  using base_t::fomSystem_;
  using base_t::fomState_;
  mutable fom_mm_action_result_type fomMMAction_;
};

//
// system implementing:
// d hat{y}/dt = MJOP f
//
template <
  class IndVarType,
  class ReducedStateType,
  class GalerkinRhsType,
  class TrialSpaceType,
  class FomSystemType,
  class HypRedOpType
  >
class GalerkinUnsteadyHypRedOdeSystem
{
public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using right_hand_side_type      = GalerkinRhsType;

  GalerkinUnsteadyHypRedOdeSystem() = delete;

  GalerkinUnsteadyHypRedOdeSystem(const TrialSpaceType & space,
				  const FomSystemType & fomSystem,
				  const HypRedOpType & hrOp)
    : space_(space),
      fomSystem_(fomSystem),
      hrOp_(hrOp),
      fomState_(fomSystem.createState()),
      fomRhs_(fomSystem.createRightHandSide())
  {}

public:
  state_type createState() const{
    return space_.get().createReducedState();
  }

  right_hand_side_type createRightHandSide() const{
    // here we assume that the action of hyOp does not
    // produce a reduced residual different than the number of basis.
    // to be precise, here we should do compute: R = MJOP f(phi x)
    const auto & phi = space_.get().viewBasis();
    return impl::CreateGalerkinRhs<right_hand_side_type>()(phi);
  }

  void rightHandSide(const state_type & galerkinState,
		     const IndVarType & timeForRhsEvaluation,
		     right_hand_side_type & galerkinRhs) const
  {
    space_.get().mapFromReducedState(galerkinState, fomState_);
    fomSystem_.get().rightHandSide(fomState_, timeForRhsEvaluation, fomRhs_);
    hrOp_(fomResidual_, galerkinRhs);
  }

protected:
  std::reference_wrapper<const TrialSpaceType> space_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  std::reference_wrapper<const HypRedOpType> hrOp_;
  mutable typename FomSystemType::state_type fomState_;
  mutable typename FomSystemType::right_hand_side_type fomRhs_;
};


//
// system class implementing:
// phi^t M phi d hat{y}/dt = MJOP f
//
template <
  class IndVarType,
  class ReducedStateType,
  class GalerkinRhsType,
  class GalerkinMassMatrixType,
  class TrialSpaceType,
  class FomSystemType,
  class HypRedOpType
  >
class GalerkinUnsteadyHypRedOdeSystemWithConstantMassMatrix
  : public GalerkinUnsteadyHypRedOdeSystem<
  IndVarType, ReducedStateType, GalerkinRhsType, TrialSpaceType, FomSystemType, HypRedOpType>
{
  using base_t = GalerkinUnsteadyHypRedOdeSystem<
    IndVarType, ReducedStateType, GalerkinRhsType,
    TrialSpaceType, FomSystemType, HypRedOpType>;

public:
  // required aliases
  using typename base_t::independent_variable_type;
  using typename base_t::state_type;
  using typename base_t::right_hand_side_type;
  using mass_matrix_type  = GalerkinMassMatrixType;

  GalerkinUnsteadyHypRedOdeSystemWithConstantMassMatrix() = delete;

  GalerkinUnsteadyHypRedOdeSystemWithConstantMassMatrix(const TrialSpaceType & space,
						const FomSystemType & fomSystem,
						const HypRedOpType & hrOp)
    : base_t(space, fomSystem, hrOp){}

public:
  mass_matrix_type createMassMatrix() const{
    const auto & phi = space_.get().viewBasis();
    return impl::CreateGalerkinMassMatrix<mass_matrix_type>()(phi);
  }

  void massMatrix(mass_matrix_type & galerkinMM) const
  {
    const auto & phi = space_.get().viewBasis();

    auto fomMMAction = fomSystem.createApplyMassMatrixResult(phi);
    fomSystem_.get().applyMassMatrix(phi, fomMMAction);

    using phi_type = typename TrialSpaceType::basis_type;
    using scalar_t = typename ::pressio::Traits<phi_type>::scalar_type;
    using cnst = ::pressio::utils::Constants<scalar_t>;
    ::pressio::ops::product(::pressio::transpose(),
			    ::pressio::nontranspose(),
			    cnst::one(), phi, fomMMAction,
			    cnst::zero(), galerkinMM);
  }

private:
  using base_t::space_;
  using base_t::fomSystem_;
  using base_t::fomState_;
};


//
// system implementing:
// d hat{y}/dt = MJOP masked(f)
//
template <
  class IndVarType,
  class ReducedStateType,
  class GalerkinRhsType,
  class TrialSpaceType,
  class FomSystemType,
  class MaskerType,
  class HypRedOpType
  >
class GalerkinUnsteadyMaskedOdeSystem
{
  using unmasked_fom_rhs_type = typename FomSystemType::right_hand_side_type;

  using masked_fom_rhs_type =
    decltype(std::declval<MaskerType const>().createApplyMaskResult
	     (std::declval<unmasked_fom_rhs_type const &>() )
	     );

public:
  // required aliases
  using independent_variable_type = IndVarType;
  using state_type                = ReducedStateType;
  using right_hand_side_type      = GalerkinRhsType;

  GalerkinUnsteadyMaskedOdeSystem() = delete;

  GalerkinUnsteadyMaskedOdeSystem(const TrialSpaceType & space,
				  const FomSystemType & fomSystem,
				  const MaskerType & masker,
				  const HypRedOpType & hrOp)
    : space_(space),
      fomSystem_(fomSystem),
      masker_(masker),
      hrOp_(hrOp),
      fomState_(fomSystem.createState()),
      unMaskedFomRhs_(fomSystem.createRightHandSide()),
      maskedFomRhs(masker.createApplyMaskResult(unMaskedFomRhs_))
  {}

public:
  state_type createState() const{
    return space_.get().createReducedState();
  }

  right_hand_side_type createRightHandSide() const{
    // here we assume that the action of hyOp does not
    // produce a reduced residual different than the number of basis.
    // to be precise, here we should do compute: MJOP masked(f)
    const auto & phi = space_.get().viewBasis();
    return impl::CreateGalerkinRhs<right_hand_side_type>()(phi);
  }

  void rightHandSide(const state_type & galerkinState,
		     const IndVarType & timeForRhsEvaluation,
		     right_hand_side_type & galerkinRhs) const
  {
    space_.get().mapFromReducedState(galerkinState, fomState_);
    fomSystem_.get().rightHandSide(fomState_, timeForRhsEvaluation, unMaskedFomRhs_);
    masker_(unMaskedFomRhs_, maskedFomRhs_);
    hrOp_(fomResidual_, galerkinRhs);
  }

protected:
  std::reference_wrapper<const TrialSpaceType> space_;
  std::reference_wrapper<const FomSystemType> fomSystem_;

  // masker and hyp-red operator
  std::reference_wrapper<const MaskerType> masker_;
  std::reference_wrapper<const HypRedOpType> hrOp_;

  mutable typename FomSystemType::state_type fomState_;
  mutable unmasked_fom_rhs_type unMaskedFomRhs_;
  mutable masked_fom_rhs_type maskedFomRhs_;
};

//
// system class implementing:
// mask(phi)^T mask(M phi) dhat{y}/dt = MJOP mask(f)
//
template <
  class IndVarType,
  class ReducedStateType,
  class GalerkinRhsType,
  class GalerkinMassMatrixType,
  class TrialSpaceType,
  class FomSystemType,
  class MaskerType,
  class HypRedOpType
  >
class GalerkinUnsteadyMaskedOdeSystemWithConstantMassMatrix
  : public GalerkinUnsteadyMaskedOdeSystem<
  IndVarType, ReducedStateType, GalerkinRhsType, TrialSpaceType,
  FomSystemType, MaskerType, HypRedOpType>
{
  using base_t = GalerkinUnsteadyMaskedOdeSystem<
    IndVarType, ReducedStateType, GalerkinRhsType,
    TrialSpaceType, FomSystemType, MaskerType, HypRedOpType>;

public:
  // required aliases
  using typename base_t::independent_variable_type;
  using typename base_t::state_type;
  using typename base_t::right_hand_side_type;
  using mass_matrix_type  = GalerkinMassMatrixType;

  GalerkinUnsteadyMaskedOdeSystemWithConstantMassMatrix() = delete;

  GalerkinUnsteadyMaskedOdeSystemWithConstantMassMatrix(const TrialSpaceType & space,
							const FomSystemType & fomSystem,
							const MaskerType & masker,
							const HypRedOpType & hrOp)
    : base_t(space, fomSystem, masker, hrOp){}

public:
  mass_matrix_type createMassMatrix() const{
    const auto & phi = space_.get().viewBasis();
    return impl::CreateGalerkinMassMatrix<mass_matrix_type>()(phi);
  }

  void massMatrix(mass_matrix_type & galerkinMM) const
  {
    const auto & phi = space_.get().viewBasis();

    // compute M*phi
    auto unMaskedFomMMAction = fomSystem.createApplyMassMatrixResult(phi);
    fomSystem_.get().applyMassMatrix(phi, unMaskedFomMMAction);

    // masked(M*phi)
    auto maskedFomMMAction = masker_.createApplyMaskResult(unMaskedFomMMAction);
    masker_(unMaskedFomMMAction, maskedFomMMAction);

    auto maskedPhi = masker_.createApplyMaskResult(phi);
    masker_(phi, maskedPhi);

    using phi_type = typename TrialSpaceType::basis_type;
    using scalar_t = typename ::pressio::Traits<phi_type>::scalar_type;
    using cnst = ::pressio::utils::Constants<scalar_t>;
    ::pressio::ops::product(::pressio::transpose(),
			    ::pressio::nontranspose(),
			    cnst::one(), maskedPhi, maskedFomMMAction,
			    cnst::zero(), galerkinMM);
  }

private:
  using base_t::space_;
  using base_t::fomSystem_;
  using base_t::fomState_;
};

}}} // end pressio::rom::impl
#endif
