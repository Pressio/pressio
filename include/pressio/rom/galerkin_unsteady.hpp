
#ifndef PRESSIO_ROM_GALERKIN_UNSTEADY_HPP_
#define PRESSIO_ROM_GALERKIN_UNSTEADY_HPP_

#include "impl/reduced_operators_helpers.hpp"
#include "impl/galerkin_unsteady_explicit_problem.hpp"
#include "impl/galerkin_unsteady_implicit_problem.hpp"
#include "impl/galerkin_unsteady_systems.hpp"

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

  // sufficient to satisfy the TrialSubspace concept since
  // the AffineSpace concept subsumes the TrialSubspace one
  static_assert(TrialSubspace<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialSubspace concept");

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

  // sufficient to satisfy the TrialSubspace concept since
  // the AffineSpace concept subsumes the TrialSubspace one
  static_assert(TrialSubspace<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialSubspace concept");

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
  impl::explicit_scheme_else_throw(schemeName, "galerkin_default_explicit");

  // sufficient to satisfy the TrialSubspace concept since
  // the AffineSpace concept subsumes the TrialSubspace one
  static_assert(TrialSubspace<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialSubspace concept");

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

  static_assert(UnsteadyExplicitGalerkinHyperReductionOperator<
		HyperReductionOperatorType, reduced_rhs_type>::value,
		"HyperReductionOperatorType does not meet the UnsteadyExplicitGalerkinHyperReductionOperator");

  // the "system" implements the math
  using galerkin_system = impl::GalerkinHyperReducedOdeSystemOnlyRhs<
    independent_variable_type, reduced_state_type, reduced_rhs_type,
    TrialSpaceType, FomSystemType, HyperReductionOperatorType>;

  using return_type = impl::GalerkinUnsteadyExplicitProblem<galerkin_system>;
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

  // sufficient to satisfy the TrialSubspace concept since
  // the AffineSpace concept subsumes the TrialSubspace one
  static_assert(TrialSubspace<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialSubspace concept");

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

  static_assert(UnsteadyExplicitGalerkinHyperReductionOperator<
		HyperReductionOperatorType, reduced_rhs_type>::value,
		"HyperReductionOperatorType does not meet the UnsteadyExplicitGalerkinHyperReductionOperator");

  // the "system" implements the math
  using galerkin_system = impl::GalerkinMaskedOdeSystemOnlyRhs<
    independent_variable_type, reduced_state_type, reduced_rhs_type,
    TrialSpaceType, FomSystemType, RhsMaskerType, HyperReductionOperatorType>;

  using return_type = impl::GalerkinUnsteadyExplicitProblem<galerkin_system>;
  return return_type(schemeName, trialSpace, fomSystem, rhsMasker, hrOp);
}

}}} // end pressio::rom::galerkin
#endif
