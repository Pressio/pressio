
#ifndef PRESSIO_ROM_LSPG_STEADY_HPP_
#define PRESSIO_ROM_LSPG_STEADY_HPP_

#include "./impl/lspg_steady_system_default.hpp"
#include "./impl/lspg_steady_system_masked.hpp"

namespace pressio{ namespace rom{ namespace lspg{

// impl-wise, default and hyp-red LSPG are the same
template<
  class TrialSubspaceType,
  class FomSystemType
  >
auto create_steady_problem(const TrialSubspaceType & trialSpace,
			   const FomSystemType & fomSystem)
{
  static_assert(PossiblyAffineTrialColumnSubspace<TrialSubspaceType>::value,
		"TrialSubspaceType does not meet the correct concept");

  static_assert(SteadyFomWithJacobianAction<
		FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value,
		"FomSystemType does not meet the SteadyFomWithJacobianAction concept");

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
  class ResidualMaskerType,
  class JacobianActionMaskerType
  >
auto create_steady_problem(TrialSubspaceType & trialSpace,
			   const FomSystemType & fomSystem,
			   const ResidualMaskerType & rMasker,
			   const JacobianActionMaskerType & jaMasker)
{
  static_assert(PossiblyAffineTrialColumnSubspace<TrialSubspaceType>::value,
		"TrialSubspaceType does not meet the correct concept");

  static_assert(SteadyFomWithJacobianAction<
		FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value,
		"FomSystemType does not meet the SteadyFomWithJacobianAction concept");

  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  // both maskers must be invariant in time
  static_assert(TimeInvariantMasker<ResidualMaskerType>::value,
		"ResidualMaskerType does not meet the TimeInvariantMasker concept");
  static_assert(TimeInvariantMasker<JacobianActionMaskerType>::value,
		"JacobianActionMaskerType does not meet the TimeInvariantMasker concept");

  // ensure the residual masker acts on the FOM residual type
  static_assert(std::is_same<
		typename ResidualMaskerType::operand_type,
		typename FomSystemType::residual_type>::value == true,
		"mismatching types of residual masker and fom residual");

  // ensure the jacobian action masker acts on the FOM jacobian action type
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<typename TrialSubspaceType::basis_matrix_type const &>()) );
  static_assert(std::is_same<
		typename JacobianActionMaskerType::operand_type,
		fom_jac_action_result_type>::value == true,
		"mismatching types of jacobian action masker and fom");

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using system_type = impl::LspgSteadyMaskedSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType,
    ResidualMaskerType, JacobianActionMaskerType>;
  return system_type(trialSpace, fomSystem, rMasker, jaMasker);
}

}}} // end pressio::rom::lspg
#endif
