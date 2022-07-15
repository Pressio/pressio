
#ifndef PRESSIO_ROM_ALL_IMPL_HPP_
#define PRESSIO_ROM_ALL_IMPL_HPP_

namespace pressio{ namespace rom{

// ========================================================
// ========================================================
namespace impl{
// ========================================================
// ========================================================

void explicit_scheme_else_throw(::pressio::ode::StepScheme name,
				       const std::string & str){
  if (!::pressio::ode::is_explicit_scheme(name)){
    throw std::runtime_error(str + " requires an explicit stepper");
  }
}

void implicit_scheme_else_throw(::pressio::ode::StepScheme name,
				const std::string & str){
  if (!::pressio::ode::is_implicit_scheme(name)){
    throw std::runtime_error(str + " requires an implicit stepper");
  }
}
// ------------------------------------------

// //
// // system class implementing:
// // d hat{y}/dt = phi^T f
// //
// template <
//   class IndVarType,
//   class ReducedStateType,
//   class GalerkinRhsType,
//   class ManifoldType,
//   class FomSystemType
//   >
// class DefaultGalerkinOdeSystemWithRhs
// {
// public:
//   // required aliases
//   using independent_variable_type = IndVarType;
//   using state_type                = ReducedStateType;
//   using right_hand_side_type      = GalerkinRhsType;

//   DefaultGalerkinOdeSystemWithRhs() = delete;

//   DefaultGalerkinOdeSystemWithRhs(ManifoldType & trialSpace,
// 				const FomSystemType & fomSystem)
//     : trialSpace_(trialSpace), fomSystem_(fomSystem),
//       fomState_(fomSystem.createState()),
//       fomRhs_(fomSystem.createRightHandSide()){}

// public:
//   state_type createState() const{
//     return trialSpace_.get().createReducedState();
//   }

//   right_hand_side_type createRightHandSide() const
//   {
//     const auto & phi = trialSpace_.get().viewManifoldJacobian();
//     return impl::CreateGalerkinRhs<right_hand_side_type>()(phi);
//   }

//   void rightHandSide(const state_type & galerkinState,
// 		     const IndVarType & timeForRhsEvaluation,
// 		     right_hand_side_type & galerkinRhs) const
//   {
//     trialSpace_.get().mapFromReducedState(galerkinState, fomState_);
//     fomSystem_.get().rightHandSide(fomState_, timeForRhsEvaluation, fomRhs_);

//     const auto & phi = trialSpace_.get().updateAndViewManifoldJacobian(galerkinState);

//     using phi_type = typename ManifoldType::manifold_jacobian_type;
//     using scalar_t = typename ::pressio::Traits<phi_type>::scalar_type;
//     using cnst = ::pressio::utils::Constants<scalar_t>;
//     ::pressio::ops::product(::pressio::transpose(), cnst::one(),
// 			    phi, fomRhs_, cnst::zero(), galerkinRhs);
//   }

// protected:
//   std::reference_wrapper<ManifoldType> trialSpace_;
//   std::reference_wrapper<const FomSystemType> fomSystem_;
//   mutable typename FomSystemType::state_type fomState_;
//   mutable typename FomSystemType::right_hand_side_type fomRhs_;
// };

// //
// // system class implementing:
// // phi^t M(y) phi d hat{y}/dt = phi^T f
// //
// template <
//   class IndVarType,
//   class ReducedStateType,
//   class GalerkinRhsType,
//   class GalerkinMassMatrixType,
//   class ManifoldType,
//   class FomSystemType
//   >
// class DefaultGalerkinOdeSystemWithRhsAndMassMatrix
//   : public DefaultGalerkinOdeSystemWithRhs<
//   IndVarType, ReducedStateType, GalerkinRhsType, ManifoldType, FomSystemType>
// {
//   using base_t = DefaultGalerkinOdeSystemWithRhs<
//     IndVarType, ReducedStateType, GalerkinRhsType,
//     ManifoldType, FomSystemType>;

//   // need to deduce the type of the action of the mass matrix
//   using manifold_jac_type = typename ManifoldType::manifold_jacobian_type;
//   using fom_mm_action_result_type =
//     decltype(std::declval<FomSystemType const>().createApplyMassMatrixResult
// 	     (std::declval<manifold_jac_type const &>())
// 	     );

// public:
//   // required aliases
//   using typename base_t::independent_variable_type;
//   using typename base_t::state_type;
//   using typename base_t::right_hand_side_type;
//   using mass_matrix_type  = GalerkinMassMatrixType;

//   DefaultGalerkinOdeSystemWithRhsAndMassMatrix() = delete;

//   DefaultGalerkinOdeSystemWithRhsAndMassMatrix(ManifoldType & trialSpace,
// 					      const FomSystemType & fomSystem)
//     : base_t(trialSpace, fomSystem),
//       fomMMAction_(fomSystem.createApplyMassMatrixResult
// 		   (trialSpace_.get().viewManifoldJacobian())
// 		   )
//   {}

// public:
//   mass_matrix_type createMassMatrix() const{
//     const auto & phi = trialSpace_.get().viewManifoldJacobian();
//     return impl::CreateGalerkinMassMatrix<mass_matrix_type>()(phi);
//   }

//   void massMatrix(const state_type & galerkinState,
// 		  const independent_variable_type & timeForEvaluation,
// 		  mass_matrix_type & galerkinMM) const
//   {

//     trialSpace_.get().mapFromReducedState(galerkinState, fomState_);

//     const auto & phi = trialSpace_.get().updateAndViewManifoldJacobian(galerkinState);
//     fomSystem_.get().applyMassMatrix(fomState_, phi, timeForEvaluation, fomMMAction_);

//     using phi_type = typename ManifoldType::manifold_jacobian_type;
//     using scalar_t = typename ::pressio::Traits<phi_type>::scalar_type;
//     using cnst = ::pressio::utils::Constants<scalar_t>;
//     ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
// 			    cnst::one(), phi, fomMMAction_, cnst::zero(), galerkinMM);
//   }

// private:
//   using base_t::trialSpace_;
//   using base_t::fomSystem_;
//   using base_t::fomState_;
//   mutable fom_mm_action_result_type fomMMAction_;
// };


// //
// // system class implementing:
// // d hat{y}/dt = phi^T f with FOM jacobian action
// //
// template <
//   class IndVarType,
//   class ReducedStateType,
//   class GalerkinRhsType,
//   class GalerkinJacobianType,
//   class ManifoldType,
//   class FomSystemType
//   >
// class DefaultGalerkinOdeSystemWithRhsAndJacobian
//   : public DefaultGalerkinOdeSystemWithRhs<
//   IndVarType, ReducedStateType, GalerkinRhsType, ManifoldType, FomSystemType>
// {
//   using base_t = DefaultGalerkinOdeSystemWithRhs<
//     IndVarType, ReducedStateType, GalerkinRhsType,
//     ManifoldType, FomSystemType>;

//   // need to deduce the type of the action of the fom jacobian
//   using manifold_jac_type = typename ManifoldType::manifold_jacobian_type;
//   using fom_jac_action_result_type =
//     decltype(std::declval<FomSystemType const>().createApplyJacobianResult
// 	     (std::declval<manifold_jac_type const &>())
// 	     );

// public:
//   // required aliases
//   using typename base_t::independent_variable_type;
//   using typename base_t::state_type;
//   using typename base_t::right_hand_side_type;
//   using jacobian_type  = GalerkinJacobianType;

//   DefaultGalerkinOdeSystemWithRhsAndJacobian() = delete;

//   DefaultGalerkinOdeSystemWithRhsAndJacobian(ManifoldType & trialSpace,
// 				const FomSystemType & fomSystem)
//     : base_t(trialSpace, fomSystem),
//       fomJacAction_(fomSystem.createApplyJacobianResult
// 		   (trialSpace_.get().viewManifoldJacobian())
// 		   )
//   {}

// public:
//   jacobian_type createJacobian() const{
//     const auto & phi = trialSpace_.get().viewManifoldJacobian();
//     return impl::CreateGalerkinJacobian<jacobian_type>()(phi);
//   }

//   void jacobian(const state_type & galerkinState,
// 		const independent_variable_type & timeForEvaluation,
// 		jacobian_type & galerkinJac) const
//   {

//     trialSpace_.get().mapFromReducedState(galerkinState, fomState_);

//     const auto & phi = trialSpace_.get().updateAndViewManifoldJacobian(galerkinState);
//     fomSystem_.get().applyJacobian(fomState_, phi, timeForEvaluation, fomJacAction_);

//     using phi_type = typename ManifoldType::manifold_jacobian_type;
//     using scalar_t = typename ::pressio::Traits<phi_type>::scalar_type;
//     using cnst = ::pressio::utils::Constants<scalar_t>;
//     ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(),
// 			    cnst::one(), phi, fomJacAction_,
// 			    cnst::zero(), galerkinJac);
//   }

// private:
//   using base_t::trialSpace_;
//   using base_t::fomSystem_;
//   using base_t::fomState_;
//   mutable fom_jac_action_result_type fomJacAction_;
// };

// //
// // explicit galerkin problem class
// //
// template <bool callableWithExtraArg, class ManifoldType, class GalSystem, class StepperType>
// class GalerkinExplicitProblem
// {
//   std::reference_wrapper<ManifoldType> trialSpace_;
//   GalSystem galSystem_;
//   StepperType stepper_;

// public:
//   using state_type = typename GalSystem::state_type;
//   using independent_variable_type  = typename GalSystem::independent_variable_type;

//   GalerkinExplicitProblem() = delete;

//   template<class FomSystemType>
//   GalerkinExplicitProblem(::pressio::ode::StepScheme schemeName,
// 			  ManifoldType & trialSpace,
// 			  const FomSystemType & fomSystem)
//     : trialSpace_(trialSpace),
//       galSystem_(trialSpace_, fomSystem),
//       stepper_( ::pressio::ode::create_explicit_stepper(schemeName, galSystem_) )
//   {}

//   template<
//     bool _callableWithExtraArg = callableWithExtraArg,
//     mpl::enable_if_t<!_callableWithExtraArg, int> = 0>
//   void operator()(state_type & state,
// 		  pressio::ode::StepStartAt<independent_variable_type> sStart,
// 		  pressio::ode::StepCount sCount,
// 		  pressio::ode::StepSize<independent_variable_type> sSize)
//   {
//     stepper_(state, sStart, sCount, sSize);
//   }

//   template<
//     class ExtraArg,
//     bool _callableWithExtraArg = callableWithExtraArg,
//     mpl::enable_if_t<_callableWithExtraArg, int> = 0
//     >
//   void operator()(state_type & state,
// 		  pressio::ode::StepStartAt<independent_variable_type> sStart,
// 		  pressio::ode::StepCount sCount,
// 		  pressio::ode::StepSize<independent_variable_type> sSize,
// 		  ExtraArg && extra)
//   {
//     stepper_(state, sStart, sCount, sSize,
// 	     std::forward<ExtraArg>(extra));
//   }
// };

// //
// // implicit galerkin problem class
// //
// template <class ManifoldType, class GalSystem, class StepperType>
// class GalerkinImplicitProblem
// {
//   std::reference_wrapper<ManifoldType> trialSpace_;
//   GalSystem galSystem_;
//   StepperType stepper_;

// public:
//   using independent_variable_type  = typename GalSystem::independent_variable_type;
//   using state_type = typename GalSystem::state_type;
//   using residual_type = typename GalSystem::right_hand_side_type;
//   using jacobian_type = typename GalSystem::jacobian_type;

//   GalerkinImplicitProblem() = delete;

//   template<class FomSystemType>
//   GalerkinImplicitProblem(::pressio::ode::StepScheme schemeName,
// 			  ManifoldType & trialSpace,
// 			  const FomSystemType & fomSystem)
//     : trialSpace_(trialSpace),
//       galSystem_(trialSpace_, fomSystem),
//       stepper_( ::pressio::ode::create_implicit_stepper(schemeName, galSystem_) )
//   {}

//   template<class ExtraArg>
//   void operator()(state_type & state,
// 		  pressio::ode::StepStartAt<independent_variable_type> sStart,
// 		  pressio::ode::StepCount sCount,
// 		  pressio::ode::StepSize<independent_variable_type> sSize,
// 		  ExtraArg && extra)
//   {
//     stepper_(state, sStart, sCount, sSize,
// 	     std::forward<ExtraArg>(extra));
//   }

//   residual_type createResidual() const{
//     return stepper_.createResidual();
//   }

//   jacobian_type createJacobian() const{
//     return stepper_.createJacobian();
//   }

//   void residual(const state_type & odeState, residual_type & R) const{
//     stepper_.residual(odeState, R);
//   }

//   void jacobian(const state_type & odeState, jacobian_type & J) const{
//     stepper_.jacobian(odeState, J);
//   }
// };


// ========================================================
// ========================================================
} // end impl
// ========================================================
// ========================================================


// template<
//   class ManifoldType, class FomSystemType,
//   mpl::enable_if_t<
//     SemiDiscreteFomWithRhs<FomSystemType>::value
//     && !SemiDiscreteFomWithRhsAndMassMatrixAction<
//       FomSystemType, typename ManifoldType::manifold_jacobian_type
//       >::value, int > = 0
//   >
// auto create_default_explicit_problem(::pressio::ode::StepScheme schemeName,
// 				     ManifoldType & trialSpace,
// 				     const FomSystemType & fomObj)
// {

//   impl::explicit_scheme_else_throw(schemeName, "galerkin_default_explicit");

//   using ind_var_type = typename FomSystemType::time_type;

//   // rom state and rhs have same type
//   using rom_state_type = typename ManifoldType::reduced_state_type;
//   using rom_rhs_type = rom_state_type;

//   using gal_system = impl::DefaultGalerkinOdeSystemWithRhs<
//     ind_var_type, rom_state_type, rom_rhs_type, ManifoldType, FomSystemType>;
//   using stepper_type =
//     decltype(::pressio::ode::create_explicit_stepper(schemeName,
// 						     std::declval<gal_system &>()
// 						     ));


//   using return_type = impl::GalerkinExplicitProblem<
//     false /*because no mm is present*/,
//     ManifoldType, gal_system, stepper_type>;
//   return return_type(schemeName, trialSpace, fomObj);
// }

// template<
//   class ManifoldType, class FomSystemType,
//   mpl::enable_if_t<
//     SemiDiscreteFomWithRhsAndMassMatrixAction<
//       FomSystemType, typename ManifoldType::manifold_jacobian_type
//       >::value, int > = 0
//   >
// auto create_default_explicit_problem(::pressio::ode::StepScheme schemeName,
// 				     ManifoldType & trialSpace,
// 				     const FomSystemType & fomObj)
// {

//   impl::explicit_scheme_else_throw(schemeName, "galerkin_default_explicit");
//   using ind_var_type = typename FomSystemType::time_type;

//   // rom state and rhs have same type
//   using rom_state_type = typename ManifoldType::reduced_state_type;
//   using rom_rhs_type = rom_state_type;

//   // figure out what is the mass matrix type from the state
//   using rom_mm_type = typename impl::galerkin_mass_matrix_type_if_state_is<rom_state_type>::type;

//   using gal_system = impl::DefaultGalerkinOdeSystemWithRhsAndMassMatrix<
//     ind_var_type, rom_state_type, rom_rhs_type, rom_mm_type, ManifoldType, FomSystemType>;
//   using stepper_type =
//     decltype(::pressio::ode::create_explicit_stepper(schemeName,
// 						     std::declval<gal_system &>()
// 						     ));

//   using return_type = impl::GalerkinExplicitProblem<
//     true /*because we have mm so its operator() needs a lin solver arg*/,
//     ManifoldType, gal_system, stepper_type>;
//   return return_type(schemeName, trialSpace, fomObj);
// }

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

//*********************************************************
} // end namespace galerkin
//*********************************************************


}} // end pressio::rom
#endif
