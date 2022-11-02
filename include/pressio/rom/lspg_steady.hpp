
#ifndef ROM_LSPG_STEADY_HPP_
#define ROM_LSPG_STEADY_HPP_

#include "./impl/lspg_steady_system_default.hpp"
#include "./impl/lspg_steady_system_masked.hpp"
#include "./rom_lspg_unsteady_hypred_updater_trilinos.hpp"

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
  class MaskerType
#if not defined PRESSIO_ENABLE_CXX20
  , mpl::enable_if_t<
      steady::ComposableIntoMaskedProblem<TrialSubspaceType, FomSystemType, MaskerType>::value,
      int> = 0
#endif
  >
#ifdef PRESSIO_ENABLE_CXX20
  requires steady::ComposableIntoMaskedProblem<
	TrialSubspaceType, FomSystemType, MaskerType>
#endif
auto create_steady_problem(TrialSubspaceType & trialSpace,
			   const FomSystemType & fomSystem,
			   const MaskerType & masker)
{

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using system_type = impl::LspgSteadyMaskedSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, MaskerType>;
  return system_type(trialSpace, fomSystem, masker);
}

namespace experimental{
// -------------------------------------------------------------
// default or hyp-red with preconditioning
// impl-wise, default and hyp-red LSPG are the same
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class PreconditionerType
#if not defined PRESSIO_ENABLE_CXX20
  , mpl::enable_if_t<
      steady::ComposableIntoDefaultOrHyperReducedProblem<
	TrialSubspaceType, FomSystemType>::value &&
      !steady::ComposableIntoMaskedProblem<
	TrialSubspaceType, FomSystemType, PreconditionerType>::value,
      int> = 0
#endif
  >
#ifdef PRESSIO_ENABLE_CXX20
  requires steady::ComposableIntoDefaultOrHyperReducedProblem<TrialSubspaceType, FomSystemType>
      && (!steady::ComposableIntoMaskedProblem<TrialSubspaceType, FomSystemType, PreconditionerType>)
#endif
auto create_steady_problem(const TrialSubspaceType & trialSpace,
			   const FomSystemType & fomSystem,
			   const PreconditionerType & preconditioner)
{

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using system_type = impl::LspgSteadyDefaultSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, PreconditionerType>;
  return system_type(trialSpace, fomSystem, preconditioner);
}

// -------------------------------------------------------------
// masked with preconditioning
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class PreconditionerType>
#ifdef PRESSIO_ENABLE_CXX20
  requires steady::ComposableIntoMaskedProblem<TrialSubspaceType, FomSystemType, MaskerType>
#endif
auto create_steady_problem(TrialSubspaceType & trialSpace,
			   const FomSystemType & fomSystem,
			   const MaskerType & masker,
			   const PreconditionerType & preconditioner)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert(steady::ComposableIntoMaskedProblem<
		TrialSubspaceType, FomSystemType, MaskerType>::value,
		"concept not satisfied");
#endif

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using system_type = impl::LspgSteadyMaskedSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, MaskerType, PreconditionerType>;
  return system_type(trialSpace, fomSystem, masker, preconditioner);
}

} // end experimental

}}} // end pressio::rom::lspg
#endif  // ROM_LSPG_STEADY_HPP_
