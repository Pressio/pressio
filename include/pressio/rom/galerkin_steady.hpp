
#ifndef PRESSIO_ROM_GALERKIN_STEADY_HPP_
#define PRESSIO_ROM_GALERKIN_STEADY_HPP_

#include "./impl/galerkin_helpers.hpp"
#include "./impl/galerkin_steady_system_default.hpp"
#include "./impl/galerkin_steady_system_hypred.hpp"
#include "./impl/galerkin_steady_system_masked.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

template<
  class TrialSubspaceType,
  class FomSystemType>
#ifdef PRESSIO_ENABLE_CXX20
  requires steady::ProjectableOnPossiblyAffineSubspace<FomSystemType, TrialSubspaceType>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSubspace,
			   const FomSystemType & fomSystem)
{
#if not defined PRESSIO_ENABLE_CXX20
  static_assert(steady::ProjectableOnPossiblyAffineSubspace<
		FomSystemType, TrialSubspaceType>::value,
		"DefaultProjectableWith not satisfied");
#endif

  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  using reduced_state_type    = typename TrialSubspaceType::reduced_state_type;
  using default_types         = SteadyGalerkinDefaultOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jac_type      = typename default_types::reduced_jacobian_type;

  using return_type = impl::GalerkinSteadyDefaultSystem<
    reduced_state_type, reduced_residual_type,
    reduced_jac_type, TrialSubspaceType, FomSystemType>;
  return return_type(trialSubspace, fomSystem);
}

template<
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType>
#ifdef PRESSIO_ENABLE_CXX20
  requires steady::HyperReduceableWith<FomSystemType, HyperReducerType, TrialSubspaceType>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSubspace,
			   const FomSystemType & fomSystem,
			   const HyperReducerType & hypReducer)
{
#if not defined PRESSIO_ENABLE_CXX20
  static_assert(steady::HyperReduceableWith<
		FomSystemType, HyperReducerType, TrialSubspaceType>::value,
		"steady::HyperReduceableWith not satisfied");
#endif

  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type >::value == true,
		"Mismatching fom states");

  using reduced_state_type    = typename TrialSubspaceType::reduced_state_type;
  using default_types         = SteadyGalerkinDefaultOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jac_type      = typename default_types::reduced_jacobian_type;

  using return_type = impl::GalerkinSteadyHypRedSystem<
    reduced_state_type, reduced_residual_type, reduced_jac_type,
    TrialSubspaceType, FomSystemType, HyperReducerType>;
  return return_type(trialSubspace, fomSystem, hypReducer);
}

template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType>
#ifdef PRESSIO_ENABLE_CXX20
requires steady::HyperReduceableAndMaskableWith<
            FomSystemType, MaskerType, HyperReducerType, TrialSubspaceType>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSubspace,
			   const FomSystemType & fomSystem,
			   const MaskerType & masker,
			   const HyperReducerType & hypReducer)
{
#if not defined PRESSIO_ENABLE_CXX20
  static_assert(steady::HyperReduceableAndMaskableWith<
		FomSystemType, MaskerType, HyperReducerType, TrialSubspaceType>::value,
		"HyperReduceableAndMaskableWith not satisfied");
#endif

  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  using reduced_state_type    = typename TrialSubspaceType::reduced_state_type;
  using default_types         = SteadyGalerkinDefaultOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jac_type      = typename default_types::reduced_jacobian_type;

  using return_type = impl::GalerkinSteadyMaskedSystem<
    reduced_state_type, reduced_residual_type, reduced_jac_type,
    TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType>;
  return return_type(trialSubspace, fomSystem, masker, hypReducer);
}

}}} // end pressio::rom::galerkin
#endif
