
#ifndef ROM_GALERKIN_UNSTEADY_IMPLICIT_HPP_
#define ROM_GALERKIN_UNSTEADY_IMPLICIT_HPP_

#include "impl/galerkin_helpers.hpp"
#include "impl/galerkin_unsteady_fom_states_manager.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_hypred_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_masked_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_fully_discrete_fom.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

// -------------------------------------------------------------
// fully-discrete
// -------------------------------------------------------------
template<
  std::size_t TotalNumberOfDesiredStates,
  class TrialSubspaceType,
  class FomSystemType>
#ifdef PRESSIO_ENABLE_CXX20
  requires FullyDiscreteSystemWithJacobianAction<
	FomSystemType, TotalNumberOfDesiredStates, typename TrialSubspaceType::basis_matrix_type>
#endif
auto create_unsteady_implicit_problem(const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert(FullyDiscreteSystemWithJacobianAction<
		FomSystemType, TotalNumberOfDesiredStates,
		typename TrialSubspaceType::basis_matrix_type>::value,
		"concept not satisfied");
#endif

  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type        = typename TrialSubspaceType::reduced_state_type;
  using default_types             = ImplicitGalerkinDefaultOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinDefaultFullyDiscreteSystem<
    TotalNumberOfDesiredStates, independent_variable_type,
    reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType>;

  galerkin_system galSystem(trialSpace, fomSystem);
  return ::pressio::ode::create_implicit_stepper<
    TotalNumberOfDesiredStates>(std::move(galSystem));
}

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
  static_assert(unsteadyimplicit::ComposableIntoDefaultProblem<
		TrialSubspaceType, FomSystemType>::value,
		"default concept not met");
#endif

  impl::valid_scheme_for_implicit_galerkin_else_throw(schemeName, "galerkin_default_implicit");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types      = ImplicitGalerkinDefaultOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinDefaultOdeSystemRhsAndJacobian<
    independent_variable_type, reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType>;

  galerkin_system galSystem(trialSpace, fomSystem);
  return ::pressio::ode::create_implicit_stepper(schemeName, std::move(galSystem));
}

// -------------------------------------------------------------
// hyper-reduced
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType>
#ifdef PRESSIO_ENABLE_CXX20
  requires unsteadyimplicit::ComposableIntoHyperReducedProblem<
	TrialSubspaceType, FomSystemType, HyperReducerType>
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

  impl::valid_scheme_for_implicit_galerkin_else_throw(schemeName, "galerkin_hypred_implicit");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types      = ImplicitGalerkinDefaultOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinHypRedOdeSystemRhsAndJacobian<
    independent_variable_type, reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType, HyperReducerType>;

  galerkin_system galSystem(trialSpace, fomSystem, hyperReducer);
  return ::pressio::ode::create_implicit_stepper(schemeName, std::move(galSystem));
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

  impl::valid_scheme_for_implicit_galerkin_else_throw(schemeName, "galerkin_masked_implicit");

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

  galerkin_system galSystem(trialSpace, fomSystem, masker, hyperReducer);
  return ::pressio::ode::create_implicit_stepper(schemeName, std::move(galSystem));
}

}}} // end pressio::rom::galerkin
#endif  // ROM_GALERKIN_UNSTEADY_IMPLICIT_HPP_
