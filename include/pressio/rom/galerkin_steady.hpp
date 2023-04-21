
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
requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
&& RealValuedSteadyFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
&& std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSpace,   /*(1)*/
			   const FomSystemType & fomSystem)
{

  using reduced_state_t    = typename TrialSubspaceType::reduced_state_type;
  using reduced_r_t = impl::steady_galerkin_default_reduced_residual_t<TrialSubspaceType>;
  using reduced_j_t = impl::steady_galerkin_default_reduced_jacobian_t<TrialSubspaceType>;
  using return_type = impl::GalerkinSteadyDefaultSystem<
    reduced_state_t, reduced_r_t, reduced_j_t,
    TrialSubspaceType, FomSystemType>;
  return return_type(trialSpace, fomSystem);
}

// ------------------------------------------------------------------------
// hyper-reduced
// ------------------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType>
#ifdef PRESSIO_ENABLE_CXX20
requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
&& RealValuedSteadyFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
&& std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSpace,   /*(2)*/
			   const FomSystemType & fomSystem,
			   const HyperReducerType & hyperReducer)
{

  using reduced_state_t = typename TrialSubspaceType::reduced_state_type;
  using reduced_r_t = impl::steady_galerkin_default_reduced_residual_t<TrialSubspaceType>;
  using reduced_j_t = impl::steady_galerkin_default_reduced_jacobian_t<TrialSubspaceType>;
  using return_type = impl::GalerkinSteadyHypRedSystem<
    reduced_state_t, reduced_r_t, reduced_j_t,
    TrialSubspaceType, FomSystemType, HyperReducerType>;
  return return_type(trialSpace, fomSystem, hyperReducer);
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
requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
&& RealValuedSteadyFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
&& std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSpace,   /*(3)*/
			   const FomSystemType & fomSystem,
			   const MaskerType & masker,
			   const HyperReducerType & hyperReducer)
{

  using reduced_state_t = typename TrialSubspaceType::reduced_state_type;
  using reduced_r_t = impl::steady_galerkin_default_reduced_residual_t<TrialSubspaceType>;
  using reduced_j_t = impl::steady_galerkin_default_reduced_jacobian_t<TrialSubspaceType>;
  using return_type = impl::GalerkinSteadyMaskedSystem<
    reduced_state_t, reduced_r_t, reduced_j_t,
    TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType>;
  return return_type(trialSpace, fomSystem, masker, hyperReducer);
}

}}} // end pressio::rom::galerkin
#endif  // ROM_GALERKIN_STEADY_HPP_
