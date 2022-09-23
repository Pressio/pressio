
#ifndef PRESSIO_ROM_LSPG_STEADY_HPP_
#define PRESSIO_ROM_LSPG_STEADY_HPP_

#include "./impl/lspg_steady_system_default.hpp"
#include "./impl/lspg_steady_system_masked.hpp"

namespace pressio{ namespace rom{ namespace lspg{

// impl-wise, default and hyp-red LSPG are the same
template<
  class TrialSubspaceType,
  class FomSystemType>
#ifdef PRESSIO_ENABLE_CXX20
requires SteadyFomWithJacobianAction<FomSystemType, TrialSubspaceType>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSpace,
			   const FomSystemType & fomSystem)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert(SteadyFomWithJacobianAction<
		FomSystemType, TrialSubspaceType>::value,
		"FomSystemType does not meet the SteadyFomWithJacobianAction concept");
#endif

  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using system_type = impl::LspgSteadyDefaultSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType>;
  return system_type(trialSpace, fomSystem);
}

template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType>
#ifdef PRESSIO_ENABLE_CXX20
requires steady::MaskableWith<FomSystemType, MaskerType, TrialSubspaceType>
#endif
auto create_steady_problem(TrialSubspaceType & trialSpace,
			   const FomSystemType & fomSystem,
			   const MaskerType & masker)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert(steady::MaskableWith<FomSystemType, MaskerType, TrialSubspaceType>::value,
		"MaskableWith not satisfied");
#endif

  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using system_type = impl::LspgSteadyMaskedSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, MaskerType>;
  return system_type(trialSpace, fomSystem, masker);
}

}}} // end pressio::rom::lspg
#endif
