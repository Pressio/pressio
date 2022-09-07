
#ifndef PRESSIO_ROM_GALERKIN_UNSTEADY_HPP_
#define PRESSIO_ROM_GALERKIN_UNSTEADY_HPP_

#include "impl/reduced_operators_helpers.hpp"
#include "impl/galerkin_unsteady_explicit_problem.hpp"
#include "impl/galerkin_unsteady_implicit_problem.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_only.hpp"
#include "impl/galerkin_unsteady_system_hypred_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_hypred_rhs_only.hpp"
#include "impl/galerkin_unsteady_system_masked_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_masked_rhs_only.hpp"

/*
  below we abuse things by using in some cases static asserts
  to check constraints even if this is not fully correct
  because constraints should impact overload resolution:
    https://timsong-cpp.github.io/cppwp/n4861/structure#footnote-154

  Since we cannot yet use c++20 concepts, we should enforce these
  constraints via e.g. SFINAE but that would yield bad error messages.
  So for now we decide to use static asserts to have readable error messages.
*/

namespace pressio{ namespace rom{

namespace impl{
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
}//end impl

namespace galerkin{

// -------------------------------------------------------------
// default
// -------------------------------------------------------------
template<
  class TrialSpaceType,
  class FomSystemType
  >
auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,
				      const TrialSpaceType & trialSpace,
				      const FomSystemType & fomSystem)
{
  impl::explicit_scheme_else_throw(schemeName, "galerkin_default_explicit");

  // sufficient to satisfy the TrialColumnSubspaceConcept concept since
  // the AffineSpace concept subsumes the TrialColumnSubspaceConcept one
  static_assert(TrialColumnSubspaceConcept<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialColumnSubspaceConcept concept");

  // sufficient to satisfy the SemiDiscreteFom concept since
  // the SemiDiscreteFomWithJacobianAction concept subsumes it
  static_assert(SemiDiscreteFom<FomSystemType>::value,
		"FomSystemType does not meet the SemiDiscreteFom concept");

  static_assert(std::is_same<typename TrialSpaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using reduced_rhs_type = reduced_state_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinDefaultOdeSystemOnlyRhs<
    independent_variable_type, reduced_state_type, reduced_rhs_type,
    TrialSpaceType, FomSystemType>;

  // a Galerkin problem contains (beside other things) a pressio stepper.
  // the reason for this is that a problem can potentially expose
  // more methods than what the underlying stepper does.
  using return_type = impl::GalerkinUnsteadyExplicitProblem<galerkin_system>;
  return return_type(schemeName, trialSpace, fomSystem);
}

template<
  class TrialSpaceType,
  class FomSystemType
  >
auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,
				      const TrialSpaceType & trialSpace,
				      const FomSystemType & fomSystem)
{
  impl::implicit_scheme_else_throw(schemeName, "galerkin_default_explicit");

  // sufficient to satisfy the TrialColumnSubspaceConcept concept since
  // the AffineSpace concept subsumes the TrialColumnSubspaceConcept one
  static_assert(TrialColumnSubspaceConcept<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialColumnSubspaceConcept concept");

  // sufficient to satisfy the SemiDiscreteFom concept since
  // the SemiDiscreteFomWithJacobianAction concept subsumes it
  static_assert(SemiDiscreteFomWithJacobianAction<
		FomSystemType, typename TrialSpaceType::basis_type>::value,
		"FomSystemType does not meet the SemiDiscreteFomWithJacobianAction concept");

  static_assert(std::is_same<typename TrialSpaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using reduced_residual_type = reduced_state_type;
  using reduced_jac_type = typename
    impl::determine_galerkin_jacobian_type_from_state<reduced_state_type>::type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinDefaultOdeSystemRhsAndJacobian<
    independent_variable_type, reduced_state_type, reduced_residual_type,
    reduced_jac_type, TrialSpaceType, FomSystemType>;

  // a Galerkin problem contains (beside other things) a pressio stepper.
  // the reason for this is that a problem can potentially expose
  // more methods than what the underlying stepper does.
  using return_type = impl::GalerkinUnsteadyImplicitProblem<galerkin_system>;
  return return_type(schemeName, trialSpace, fomSystem);
}


// -------------------------------------------------------------
// hyper-reduced
// -------------------------------------------------------------
template<
  class TrialSpaceType,
  class FomSystemType,
  class HyperReductionOperatorType
  >
auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,
				      const TrialSpaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const HyperReductionOperatorType & hrOp)
{
  impl::explicit_scheme_else_throw(schemeName, "galerkin_default_implicit");

  // sufficient to satisfy the TrialColumnSubspaceConcept concept since
  // the AffineSpace concept subsumes the TrialColumnSubspaceConcept one
  static_assert(TrialColumnSubspaceConcept<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialColumnSubspaceConcept concept");

  // sufficient to satisfy the SemiDiscreteFom concept since
  // the SemiDiscreteFomWithJacobianAction concept subsumes it
  static_assert(SemiDiscreteFom<FomSystemType>::value,
		"FomSystemType does not meet the SemiDiscreteFom concept");

  static_assert(std::is_same<typename TrialSpaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using reduced_rhs_type = reduced_state_type;

  static_assert(ExplicitGalerkinHyperReducer<
		HyperReductionOperatorType, reduced_rhs_type>::value,
		"HyperReductionOperatorType does not meet the ExplicitGalerkinHyperReducer");

  // the "system" implements the math
  using galerkin_system = impl::GalerkinHyperReducedOdeSystemOnlyRhs<
    independent_variable_type, reduced_state_type, reduced_rhs_type,
    TrialSpaceType, FomSystemType, HyperReductionOperatorType>;

  using return_type = impl::GalerkinUnsteadyExplicitProblem<galerkin_system>;
  return return_type(schemeName, trialSpace, fomSystem, hrOp);
}

template<
  class TrialSpaceType,
  class FomSystemType,
  class HyperReductionOperatorType
  >
auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,
				      const TrialSpaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const HyperReductionOperatorType & hrOp)
{
  impl::implicit_scheme_else_throw(schemeName, "galerkin_hypred_implicit");

  // sufficient to satisfy the TrialColumnSubspaceConcept concept since
  // the AffineSpace concept subsumes the TrialColumnSubspaceConcept one
  static_assert(TrialColumnSubspaceConcept<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialColumnSubspaceConcept concept");

  // sufficient to satisfy the SemiDiscreteFom concept since
  // the SemiDiscreteFomWithJacobianAction concept subsumes it
  static_assert(SemiDiscreteFomWithJacobianAction<
		FomSystemType, typename TrialSpaceType::basis_type>::value,
		"FomSystemType does not meet the SemiDiscreteFomWithJacobianAction concept");

  static_assert(std::is_same<typename TrialSpaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using reduced_residual_type = reduced_state_type;
  using reduced_jac_type = typename
    impl::determine_galerkin_jacobian_type_from_state<reduced_state_type>::type;

  static_assert(ImplicitGalerkinHyperReducer<
		HyperReductionOperatorType, reduced_residual_type, reduced_jac_type>::value,
		"HyperReductionOperatorType does not meet the ImplicitGalerkinHyperReducer");

  // the "system" implements the math
  using galerkin_system = impl::GalerkinHypRedOdeSystemRhsAndJacobian<
    independent_variable_type, reduced_state_type, reduced_residual_type,
    reduced_jac_type, TrialSpaceType, FomSystemType, HyperReductionOperatorType>;

  // a Galerkin problem contains (beside other things) a pressio stepper.
  // the reason for this is that a problem can potentially expose
  // more methods than what the underlying stepper does.
  using return_type = impl::GalerkinUnsteadyImplicitProblem<galerkin_system>;
  return return_type(schemeName, trialSpace, fomSystem, hrOp);
}

// -------------------------------------------------------------
// masked
// -------------------------------------------------------------
template<
  class TrialSpaceType,
  class FomSystemType,
  class RhsMaskerType,
  class HyperReductionOperatorType
  >
auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,
				      const TrialSpaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const RhsMaskerType & rhsMasker,
				      const HyperReductionOperatorType & hrOp)
{
  impl::explicit_scheme_else_throw(schemeName, "galerkin_default_explicit");

  // sufficient to satisfy the TrialColumnSubspaceConcept concept since
  // the AffineSpace concept subsumes the TrialColumnSubspaceConcept one
  static_assert(TrialColumnSubspaceConcept<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialColumnSubspaceConcept concept");

  // sufficient to satisfy the SemiDiscreteFom concept since
  // the SemiDiscreteFomWithJacobianAction concept subsumes it
  static_assert(SemiDiscreteFom<FomSystemType>::value,
		"FomSystemType does not meet the SemiDiscreteFom concept");

  static_assert(std::is_same<typename TrialSpaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  // masker must be invariant in time
  static_assert(TimeInvariantMasker<RhsMaskerType>::value,
		"RhsMaskerType does not meet the TimeInvariantMasker concept");

  // ensure the masker acts on the FOM rhs type
  static_assert(std::is_same<
		typename RhsMaskerType::operand_type,
		typename FomSystemType::right_hand_side_type>::value == true,
		"mismatching types of rhs masker and fom right_hand_side_type");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using reduced_rhs_type = reduced_state_type;

  static_assert(ExplicitGalerkinHyperReducer<
		HyperReductionOperatorType, reduced_rhs_type>::value,
		"HyperReductionOperatorType does not meet the ExplicitGalerkinHyperReducer");

  // the "system" implements the math
  using galerkin_system = impl::GalerkinMaskedOdeSystemOnlyRhs<
    independent_variable_type, reduced_state_type, reduced_rhs_type,
    TrialSpaceType, FomSystemType, RhsMaskerType, HyperReductionOperatorType>;

  using return_type = impl::GalerkinUnsteadyExplicitProblem<galerkin_system>;
  return return_type(schemeName, trialSpace, fomSystem, rhsMasker, hrOp);
}

template<
  class TrialSpaceType,
  class FomSystemType,
  class RhsMaskerType,
  class JacobianActionMaskerType,
  class HyperReductionOperatorType
  >
auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,
				      const TrialSpaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const RhsMaskerType & rhsMasker,
				      const JacobianActionMaskerType & jaMasker,
				      const HyperReductionOperatorType & hrOp)
{
  impl::implicit_scheme_else_throw(schemeName, "galerkin_hypred_implicit");

  // sufficient to satisfy the TrialColumnSubspaceConcept concept since
  // the AffineSpace concept subsumes the TrialColumnSubspaceConcept one
  static_assert(TrialColumnSubspaceConcept<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialColumnSubspaceConcept concept");

  // sufficient to satisfy the SemiDiscreteFom concept since
  // the SemiDiscreteFomWithJacobianAction concept subsumes it
  static_assert(SemiDiscreteFomWithJacobianAction<
		FomSystemType, typename TrialSpaceType::basis_type>::value,
		"FomSystemType does not meet the SemiDiscreteFomWithJacobianAction concept");

  static_assert(std::is_same<typename TrialSpaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  // both maskers must be invariant in time
  static_assert(TimeInvariantMasker<RhsMaskerType>::value,
		"RhsMaskerType does not meet the TimeInvariantMasker concept");
  static_assert(TimeInvariantMasker<JacobianActionMaskerType>::value,
		"JacobianActionMaskerType does not meet the TimeInvariantMasker concept");

  // ensure the rhs masker acts on the FOM rhs type
  static_assert(std::is_same<
		typename RhsMaskerType::operand_type,
		typename FomSystemType::right_hand_side_type>::value == true,
		"mismatching types of residual masker and fom rhs");

  // ensure the jacobian action masker acts on the FOM jacobian action type
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<typename TrialSpaceType::basis_type const &>()) );
  static_assert(std::is_same<
		typename JacobianActionMaskerType::operand_type,
		fom_jac_action_result_type>::value == true,
		"mismatching types of jacobian action masker and fom residual");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using reduced_residual_type = reduced_state_type;
  using reduced_jac_type = typename
    impl::determine_galerkin_jacobian_type_from_state<reduced_state_type>::type;

  static_assert(ImplicitGalerkinHyperReducer<
		HyperReductionOperatorType, reduced_residual_type, reduced_jac_type>::value,
		"HyperReductionOperatorType does not meet the ImplicitGalerkinHyperReducer");

  // the hr operator must operate on types that are consistent with
  // those deducible from the masker
  static_assert(std::is_same<
		typename HyperReductionOperatorType::right_hand_side_operand_type,
		typename RhsMaskerType::result_type>::value == true,
		"mismatching types");
  static_assert(std::is_same<
		typename HyperReductionOperatorType::jacobian_action_operand_type,
		typename JacobianActionMaskerType::result_type>::value == true,
		"mismatching types");

  // the "system" implements the math
  using galerkin_system = impl::GalerkinMaskedOdeSystemRhsAndJacobian<
    independent_variable_type, reduced_state_type, reduced_residual_type,
    reduced_jac_type, TrialSpaceType, FomSystemType,
    RhsMaskerType, JacobianActionMaskerType, HyperReductionOperatorType>;

  // a Galerkin problem contains (beside other things) a pressio stepper.
  // the reason for this is that a problem can potentially expose
  // more methods than what the underlying stepper does.
  using return_type = impl::GalerkinUnsteadyImplicitProblem<galerkin_system>;
  return return_type(schemeName, trialSpace, fomSystem,
		     rhsMasker, jaMasker, hrOp);
}

}}} // end pressio::rom::galerkin
#endif
