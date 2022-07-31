
//
// system class implementing:
// d hat{y}/dt = phi^T f with FOM jacobian action
//
template <
  class IndVarType,
  class ReducedStateType,
  class GalerkinRhsType,
  class GalerkinJacobianType,
  class ManifoldType,
  class FomSystemType
  >
class DefaultGalerkinOdeSystemWithRhsAndJacobian
  : public DefaultGalerkinOdeSystemWithRhs<
  IndVarType, ReducedStateType, GalerkinRhsType, ManifoldType, FomSystemType>
{
  using base_t = DefaultGalerkinOdeSystemWithRhs<
    IndVarType, ReducedStateType, GalerkinRhsType,
    ManifoldType, FomSystemType>;

  // need to deduce the type of the action of the fom jacobian
  using manifold_jac_type = typename ManifoldType::manifold_jacobian_type;
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<manifold_jac_type const &>())
	     );

public:
  // required aliases
  using typename base_t::independent_variable_type;
  using typename base_t::state_type;
  using typename base_t::right_hand_side_type;
  using jacobian_type  = GalerkinJacobianType;

  DefaultGalerkinOdeSystemWithRhsAndJacobian() = delete;

  DefaultGalerkinOdeSystemWithRhsAndJacobian(ManifoldType & trialSpace,
				const FomSystemType & fomSystem)
    : base_t(trialSpace, fomSystem),
      fomJacAction_(fomSystem.createApplyJacobianResult
		   (trialSpace_.get().viewManifoldJacobian())
		   )
  {}

public:
  jacobian_type createJacobian() const{
    const auto & phi = trialSpace_.get().viewManifoldJacobian();
    return impl::CreateGalerkinJacobian<jacobian_type>()(phi);
  }

  void jacobian(const state_type & galerkinState,
		const independent_variable_type & timeForEvaluation,
		jacobian_type & galerkinJac) const
  {

    trialSpace_.get().mapFromReducedState(galerkinState, fomState_);

    const auto & phi = trialSpace_.get().updateAndViewManifoldJacobian(galerkinState);
    fomSystem_.get().applyJacobian(fomState_, phi, timeForEvaluation, fomJacAction_);

    using phi_type = typename ManifoldType::manifold_jacobian_type;
    using scalar_t = typename ::pressio::Traits<phi_type>::scalar_type;
    using cnst = ::pressio::utils::Constants<scalar_t>;
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
			    cnst::one(), phi, fomJacAction_,
			    cnst::zero(), galerkinJac);
  }

private:
  using base_t::trialSpace_;
  using base_t::fomSystem_;
  using base_t::fomState_;
  mutable fom_jac_action_result_type fomJacAction_;
};

// template<
//   class ManifoldType, class FomSystemType,
//   mpl::enable_if_t<
//     SemiDiscreteFomWithRhsAndJacobianAction<
//       FomSystemType, typename ManifoldType::manifold_jacobian_type>::value, int > = 0
//   >
// auto create_default_implicit_problem(::pressio::ode::StepScheme schemeName,
// 				     ManifoldType & trialSpace,
// 				     const FomSystemType & fomObj)
// {

//   impl::implicit_scheme_else_throw(schemeName, "galerkin_default_implicit");

//   using ind_var_type = typename FomSystemType::time_type;

//   // rom state and rhs have same type
//   using rom_state_type = typename ManifoldType::reduced_state_type;
//   using rom_rhs_type = rom_state_type;

//   // figure out what is the reduced jacobian type from the state
//   using rom_jac_type = typename impl::galerkin_jacobian_type_if_state_is<rom_state_type>::type;

//   using gal_system = impl::DefaultGalerkinOdeSystemWithRhsAndJacobian<
//     ind_var_type, rom_state_type, rom_rhs_type, rom_jac_type, ManifoldType, FomSystemType>;
//   using stepper_type =
//     decltype(::pressio::ode::create_implicit_stepper(schemeName,
// 						     std::declval<gal_system &>()
// 						     ));

//   using return_type = impl::GalerkinImplicitProblem<
//     ManifoldType, gal_system, stepper_type>;
//   return return_type(schemeName, trialSpace, fomObj);
// }
