
#ifndef ROM_LSPG_STEADY_HPP_
#define ROM_LSPG_STEADY_HPP_

#include "./impl/lspg_steady_system_default.hpp"
#include "./impl/lspg_steady_system_masked.hpp"
#include "./impl/noop.hpp"

namespace pressio{ namespace rom{ namespace lspg{

// -------------------------------------------------------------
// default or hyp-red
// impl-wise, default and hyp-red LSPG are the same
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType>
#ifdef PRESSIO_ENABLE_CXX20
requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
&& RealValuedSteadyFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
&& std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSpace,   /*(1)*/
			   const FomSystemType & fomSystem)

{

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = impl::NoOperation<void>;
  using system_type = impl::LspgSteadyDefaultSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, scaler_type>;
  return system_type(trialSpace, fomSystem, scaler_type{});
}

// -------------------------------------------------------------
// masked
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType>
#ifdef PRESSIO_ENABLE_CXX20
requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
&& RealValuedSteadyFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
&& std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSpace,  /*(2)*/
			   const FomSystemType & fomSystem,
			   const MaskerType & masker)
{

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = impl::NoOperation<void>;
  using system_type = impl::LspgSteadyMaskedSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, MaskerType, scaler_type>;
  return system_type(trialSpace, fomSystem, masker, scaler_type{});
}


namespace experimental{
// -------------------------------------------------------------
// default or hyp-red with scaling
// impl-wise, default and hyp-red LSPG are the same
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class ScalingOperatorType>
#ifdef PRESSIO_ENABLE_CXX20
requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
&& RealValuedSteadyFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
&& std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSpace,  /*(3)*/
			   const FomSystemType & fomSystem,
			   const ScalingOperatorType & scaler)
{

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = std::reference_wrapper<const ScalingOperatorType>;
  using system_type = impl::LspgSteadyDefaultSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, scaler_type>;
  return system_type(trialSpace, fomSystem, scaler);
}

// -------------------------------------------------------------
// masked with scaling
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class ScalingOperatorType>
#ifdef PRESSIO_ENABLE_CXX20
requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
&& RealValuedSteadyFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
&& std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSpace,  /*(4)*/
			   const FomSystemType & fomSystem,
			   const MaskerType & masker,
			   const ScalingOperatorType & scaler)
{

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = std::reference_wrapper<const ScalingOperatorType>;
  using system_type = impl::LspgSteadyMaskedSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, MaskerType, scaler_type>;
  return system_type(trialSpace, fomSystem, masker, scaler);
}

} // end experimental

}}} // end pressio::rom::lspg
#endif  // ROM_LSPG_STEADY_HPP_
