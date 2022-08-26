
#ifndef PRESSIO_ROM_LSPG_STEADY_HPP_
#define PRESSIO_ROM_LSPG_STEADY_HPP_

#include "./impl/lspg_steady_default_system.hpp"
#include "./impl/lspg_steady_masked_system.hpp"

namespace pressio{ namespace rom{ namespace lspg{

// impl-wise, default and hyp-red LSPG are the same
template<
  class TrialSpaceType,
  class FomSystemType
  >
auto create_steady_problem(const TrialSpaceType & trialSpace,
			   const FomSystemType & fomSystem)
{

  // sufficient to satisfy the TrialSubspace concept since
  // the AffineSpace concept subsumes the TrialSubspace one
  static_assert(TrialSubspace<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialSubspace concept");

  static_assert(SteadyFomWithJacobianAction<
		FomSystemType, typename TrialSpaceType::basis_type>::value,
		"FomSystemType does not meet the SteadyFomWithJacobianAction concept");

  static_assert(std::is_same<typename TrialSpaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using system_type = impl::LspgSteadyDefaultSystem<
    reduced_state_type, TrialSpaceType, FomSystemType>;
  return system_type(trialSpace, fomSystem);
}

template<
  class TrialSpaceType,
  class FomSystemType,
  class ResidualMaskerType,
  class JacobianActionMaskerType
  >
auto create_steady_problem(TrialSpaceType & trialSpace,
			   const FomSystemType & fomSystem,
			   const ResidualMaskerType & rMasker,
			   const JacobianActionMaskerType & jaMasker)
{

  // sufficient to satisfy the TrialSubspace concept since
  // the AffineSpace concept subsumes the TrialSubspace one
  static_assert(TrialSubspace<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialSubspace concept");

  static_assert(SteadyFomWithJacobianAction<
		FomSystemType, typename TrialSpaceType::basis_type>::value,
		"FomSystemType does not meet the SteadyFomWithJacobianAction concept");

  static_assert(std::is_same<typename TrialSpaceType::full_state_type,
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
	     (std::declval<typename TrialSpaceType::basis_type const &>()) );
  static_assert(std::is_same<
		typename JacobianActionMaskerType::operand_type,
		fom_jac_action_result_type>::value == true,
		"mismatching types of jacobian action masker and fom residual");

  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using system_type = impl::LspgSteadyMaskedSystem<
    reduced_state_type, TrialSpaceType, FomSystemType,
    ResidualMaskerType, JacobianActionMaskerType>;
  return system_type(trialSpace, fomSystem, rMasker, jaMasker);
}

// //
// // overload set for preconditioned problems
// // (might not need these if the preconditioner is instead provided to solver)
// //
// template<
//   class TrialSpaceType, class FomSystemType, class PreconditionerType,
//   mpl::enable_if_t<
//     SteadyFomWithJacobianAction<
//       FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
//   >
// auto create_preconditioned_problem(TrialSpaceType & trialSpace,
//            const FomSystemType & fomSystem,
//            const PreconditionerType & preconditioner)
// {
// using reduced_state_type = typename TrialSpaceType::reduced_state_type;
// using system_type = impl::LspgSteadyPreconditionedSystem<
//   reduced_state_type, TrialSpaceType, FomSystemType, PreconditionerType>;
// return system_type(trialSpace, fomSystem, preconditioner);
// }

}}} // end pressio::rom::galerkin
#endif
