
#ifndef ROM_GALERKIN_STEADY_HPP_
#define ROM_GALERKIN_STEADY_HPP_

#include "./impl/galerkin_helpers.hpp"
#include "./impl/galerkin_steady_system_default.hpp"
#include "./impl/galerkin_steady_system_hypred.hpp"
#include "./impl/galerkin_steady_system_masked.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

// ------------------------------------------------------------------------
// default
// ------------------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType>
#ifdef PRESSIO_ENABLE_CXX20
  requires steady::ComposableIntoDefaultProblem<TrialSubspaceType, FomSystemType>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSubspace,
			   const FomSystemType & fomSystem)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert(steady::ComposableIntoDefaultProblem<
		TrialSubspaceType, FomSystemType>::value,
		"steady::ComposableIntoDefaultProblem not satisfied");
#endif

  using reduced_state_type    = typename TrialSubspaceType::reduced_state_type;
  using default_types         = SteadyGalerkinDefaultOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jac_type      = typename default_types::reduced_jacobian_type;

  using return_type = impl::GalerkinSteadyDefaultSystem<
    reduced_state_type, reduced_residual_type, reduced_jac_type,
    TrialSubspaceType, FomSystemType>;
  return return_type(trialSubspace, fomSystem);
}

// ------------------------------------------------------------------------
// hyper-reduced
// ------------------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType>
#ifdef PRESSIO_ENABLE_CXX20
  requires steady::ComposableIntoHyperReducedProblem<
	TrialSubspaceType, FomSystemType, HyperReducerType>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSubspace,
			   const FomSystemType & fomSystem,
			   const HyperReducerType & hyperReducer)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert(steady::ComposableIntoHyperReducedProblem<
		TrialSubspaceType, FomSystemType, HyperReducerType>::value,
		"steady::ComposableIntoHyperReducedProblem not satisfied");
#endif

  using reduced_state_type    = typename TrialSubspaceType::reduced_state_type;
  using default_types         = SteadyGalerkinDefaultOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jac_type      = typename default_types::reduced_jacobian_type;

  using return_type = impl::GalerkinSteadyHypRedSystem<
    reduced_state_type, reduced_residual_type, reduced_jac_type,
    TrialSubspaceType, FomSystemType, HyperReducerType>;
  return return_type(trialSubspace, fomSystem, hyperReducer);
}

// ------------------------------------------------------------------------
// hyper-reduced masked
// ------------------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType>
#ifdef PRESSIO_ENABLE_CXX20
  requires steady::ComposableIntoHyperReducedMaskedProblem<
	TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSubspace,
			   const FomSystemType & fomSystem,
			   const MaskerType & masker,
			   const HyperReducerType & hyperReducer)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert(steady::ComposableIntoHyperReducedMaskedProblem<
		TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType>::value,
		"steady::ComposableIntoHyperReducedMaskedProblem not satisfied");
#endif

  using reduced_state_type    = typename TrialSubspaceType::reduced_state_type;
  using default_types         = SteadyGalerkinDefaultOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jac_type      = typename default_types::reduced_jacobian_type;

  using return_type = impl::GalerkinSteadyMaskedSystem<
    reduced_state_type, reduced_residual_type, reduced_jac_type,
    TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType>;
  return return_type(trialSubspace, fomSystem, masker, hyperReducer);
}

}}} // end pressio::rom::galerkin
#endif  // ROM_GALERKIN_STEADY_HPP_
