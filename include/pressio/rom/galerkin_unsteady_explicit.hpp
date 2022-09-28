
#ifndef PRESSIO_ROM_GALERKIN_UNSTEADY_EXPLICIT_HPP_
#define PRESSIO_ROM_GALERKIN_UNSTEADY_EXPLICIT_HPP_

#include "impl/galerkin_helpers.hpp"
#include "impl/galerkin_unsteady_explicit_problem.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_only.hpp"
#include "impl/galerkin_unsteady_system_hypred_rhs_only.hpp"
#include "impl/galerkin_unsteady_system_masked_rhs_only.hpp"

namespace pressio{ namespace rom{

namespace impl{
void explicit_scheme_else_throw(::pressio::ode::StepScheme name,
				const std::string & str){
  if (!::pressio::ode::is_explicit_scheme(name)){
    throw std::runtime_error(str + " requires an explicit stepper");
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
requires unsteadyexplicit::ComposableIntoDefaultProblem<TrialSubspaceType, FomSystemType>
#endif
auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem)
{
#if not defined PRESSIO_ENABLE_CXX20
  static_assert(unsteadyexplicit::ComposableIntoDefaultProblem<TrialSubspaceType, FomSystemType>::value,
		"default concept not met");
#endif

  impl::explicit_scheme_else_throw(schemeName, "galerkin_default_explicit");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types = ExplicitGalerkinDefaultOperatorsTraits<reduced_state_type>;
  using reduced_rhs_type = typename default_types::reduced_right_hand_side_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinDefaultOdeSystemOnlyRhs<
    independent_variable_type, reduced_state_type, reduced_rhs_type,
    TrialSubspaceType, FomSystemType>;

  // a Galerkin problem contains (beside other things) a pressio stepper.
  // the reason for this is that a problem can potentially expose
  // more methods than what the underlying stepper does.
  using return_type = impl::GalerkinUnsteadyExplicitProblem<galerkin_system>;
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
requires unsteadyexplicit::ComposableIntoHyperReducedProblem<TrialSubspaceType, FomSystemType, HyperReducerType>
#endif
auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const HyperReducerType & hyperReducer)
{
#if not defined PRESSIO_ENABLE_CXX20
  static_assert(unsteadyexplicit::ComposableIntoHyperReducedProblem<
		TrialSubspaceType, FomSystemType, HyperReducerType>::value,
		"concept not satisfied");
#endif

  impl::explicit_scheme_else_throw(schemeName, "galerkin_hypred_explicit");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types = ExplicitGalerkinDefaultOperatorsTraits<reduced_state_type>;
  using reduced_rhs_type = typename default_types::reduced_right_hand_side_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinHyperReducedOdeSystemOnlyRhs<
    independent_variable_type, reduced_state_type, reduced_rhs_type,
    TrialSubspaceType, FomSystemType, HyperReducerType>;

  using return_type = impl::GalerkinUnsteadyExplicitProblem<galerkin_system>;
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
requires unsteadyexplicit::ComposableIntoHyperReducedMaskedProblem<
  TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType>
#endif
auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const MaskerType & masker,
				      const HyperReducerType & hyperReducer)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert(unsteadyexplicit::ComposableIntoHyperReducedMaskedProblem<
		TrialSubspaceType, FomSystemType,
		MaskerType, HyperReducerType>::value,
		"masked concept not satisfied");
#endif

  impl::explicit_scheme_else_throw(schemeName, "galerkin_masked_explicit");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using reduced_rhs_type = reduced_state_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinMaskedOdeSystemOnlyRhs<
    independent_variable_type, reduced_state_type, reduced_rhs_type,
    TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType>;

  using return_type = impl::GalerkinUnsteadyExplicitProblem<galerkin_system>;
  return return_type(schemeName, trialSpace, fomSystem, masker, hyperReducer);
}

}}} // end pressio::rom::galerkin
#endif
