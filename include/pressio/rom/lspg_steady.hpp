
#ifndef ROM_LSPG_STEADY_HPP_
#define ROM_LSPG_STEADY_HPP_

#include "./impl/lspg_steady_system_default.hpp"
#include "./impl/lspg_steady_system_masked.hpp"

namespace pressio{ namespace rom{ namespace lspg{

// -------------------------------------------------------------
// default or hyp-red
// impl-wise, default and hyp-red LSPG are the same
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType>
#ifdef PRESSIO_ENABLE_CXX20
  requires steady::ComposableIntoDefaultOrHyperReducedProblem<TrialSubspaceType, FomSystemType>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSpace,
			   const FomSystemType & fomSystem)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert(steady::ComposableIntoDefaultOrHyperReducedProblem<
		TrialSubspaceType, FomSystemType>::value,
		"concept not met");
#endif

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using system_type = impl::LspgSteadyDefaultSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType>;
  return system_type(trialSpace, fomSystem);
}

// -------------------------------------------------------------
// masked
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType>
#ifdef PRESSIO_ENABLE_CXX20
  requires steady::ComposableIntoMaskedProblem<TrialSubspaceType, FomSystemType, MaskerType>
#endif
auto create_steady_problem(TrialSubspaceType & trialSpace,
			   const FomSystemType & fomSystem,
			   const MaskerType & masker)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert(steady::ComposableIntoMaskedProblem<
		TrialSubspaceType, FomSystemType, MaskerType>::value,
		"concept not satisfied");
#endif

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using system_type = impl::LspgSteadyMaskedSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, MaskerType>;
  return system_type(trialSpace, fomSystem, masker);
}

}}} // end pressio::rom::lspg
#endif  // ROM_LSPG_STEADY_HPP_
