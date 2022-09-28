
#ifndef PRESSIO_ROM_GALERKIN_UNSTEADY_IMPLICIT_HPP_
#define PRESSIO_ROM_GALERKIN_UNSTEADY_IMPLICIT_HPP_

#include "impl/galerkin_helpers.hpp"
#include "impl/galerkin_unsteady_implicit_problem.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_hypred_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_masked_rhs_and_jacobian.hpp"

namespace pressio{ namespace rom{

namespace impl{
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
  class TrialSubspaceType,
  class FomSystemType>
#ifdef PRESSIO_ENABLE_CXX20
requires unsteadyimplicit::ComposableIntoDefaultProblem<TrialSubspaceType, FomSystemType>
#endif
auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem)
{
#if not defined PRESSIO_ENABLE_CXX20
  static_assert(unsteadyimplicit::ComposableIntoDefaultProblem<TrialSubspaceType, FomSystemType>::value,
		"default concept not met");
#endif

  impl::implicit_scheme_else_throw(schemeName, "galerkin_default_implicit");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types      = ImplicitGalerkinDefaultOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinDefaultOdeSystemRhsAndJacobian<
    independent_variable_type, reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType>;

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
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType>
#ifdef PRESSIO_ENABLE_CXX20
requires unsteadyimplicit::ComposableIntoHyperReducedProblem<TrialSubspaceType, FomSystemType, HyperReducerType>
#endif
auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const HyperReducerType & hyperReducer)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert(unsteadyimplicit::ComposableIntoHyperReducedProblem<
		TrialSubspaceType, FomSystemType, HyperReducerType>::value,
		"concept not satisfied");
#endif

  impl::implicit_scheme_else_throw(schemeName, "galerkin_hypred_implicit");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types      = ImplicitGalerkinDefaultOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinHypRedOdeSystemRhsAndJacobian<
    independent_variable_type, reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType, HyperReducerType>;

  // a Galerkin problem contains (beside other things) a pressio stepper.
  // the reason for this is that a problem can potentially expose
  // more methods than what the underlying stepper does.
  using return_type = impl::GalerkinUnsteadyImplicitProblem<galerkin_system>;
  return return_type(schemeName, trialSpace, fomSystem, hyperReducer);
}

// -------------------------------------------------------------
// masked
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType>
#ifdef PRESSIO_ENABLE_CXX20
requires unsteadyimplicit::ComposableIntoHyperReducedMaskedProblem<
  TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType>
#endif
auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const MaskerType & masker,
				      const HyperReducerType & hyperReducer)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert(unsteadyimplicit::ComposableIntoHyperReducedMaskedProblem<
		TrialSubspaceType, FomSystemType,
		MaskerType, HyperReducerType>::value,
		"masked concept not satisfied");
#endif

  impl::implicit_scheme_else_throw(schemeName, "galerkin_masked_implicit");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types      = ImplicitGalerkinDefaultOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinMaskedOdeSystemRhsAndJacobian<
    independent_variable_type, reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType,
    MaskerType, HyperReducerType>;

  // a Galerkin problem contains (beside other things) a pressio stepper.
  // the reason for this is that a problem can potentially expose
  // more methods than what the underlying stepper does.
  using return_type = impl::GalerkinUnsteadyImplicitProblem<galerkin_system>;
  return return_type(schemeName, trialSpace, fomSystem,
		     masker, hyperReducer);
}

}}} // end pressio::rom::galerkin
#endif
