
#ifndef PRESSIO_ROM_GALERKIN_STEADY_HPP_
#define PRESSIO_ROM_GALERKIN_STEADY_HPP_

#include "impl/galerkin_steady_default_system.hpp"
#include "impl/galerkin_steady_hypred_system.hpp"
#include "impl/galerkin_steady_masked_system.hpp"

/*
  below we abuse things by using in some cases static asserts
  to check constraints even if this is not fully correct
  because constraints should impact overload resolution:
    https://timsong-cpp.github.io/cppwp/n4861/structure#footnote-154

  Since we cannot yet use c++20 concepts, we should enforce these
  constraints via e.g. SFINAE but that would yield bad error messages.
  So for now we decide to use static asserts to have readable error messages.
*/

namespace pressio{ namespace rom{ namespace galerkin{

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
  using reduced_residual_type = reduced_state_type;

  // figure out what is the reduced jacobian type from the state
  // for example, if the reduced state is Eigen vector,
  // it makes sense to use an Eigen dense matrix to store
  // the Galerkin jacobian since all reduced operators are dense
  using reduced_jac_type = typename
    impl::determine_galerkin_jacobian_type_from_state<reduced_state_type>::type;

  using return_type = impl::GalerkinSteadyDefaultSystem<
    reduced_state_type, reduced_residual_type,
    reduced_jac_type, TrialSpaceType, FomSystemType>;
  return return_type(trialSpace, fomSystem);
}

template<
  class TrialSpaceType,
  class FomSystemType,
  class HyperReductionOperatorType
  >
auto create_steady_problem(const TrialSpaceType & trialSpace,
			   const FomSystemType & fomSystem,
			   const HyperReductionOperatorType & hrOp)
{
  static_assert(TrialSubspace<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialSubspace concept");

  static_assert(SteadyFomWithJacobianAction<
		FomSystemType, typename TrialSpaceType::basis_type>::value,
		"FomSystemType does not meet the SteadyFomWithJacobianAction concept");

  static_assert(std::is_same<typename TrialSpaceType::full_state_type,
		typename FomSystemType::state_type >::value == true,
		"Mismatching fom states");

  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using reduced_residual_type = reduced_state_type;
  using reduced_jac_type = typename
    impl::determine_galerkin_jacobian_type_from_state<reduced_state_type>::type;

  static_assert(SteadyGalerkinHyperReducer<
		HyperReductionOperatorType, reduced_residual_type, reduced_jac_type>::value,
		"HyperReductionOperatorType does not meet the SteadyGalerkinHyperReducer");

  using return_type = impl::GalerkinSteadyHypRedSystem<
    reduced_state_type, reduced_residual_type, reduced_jac_type,
    TrialSpaceType, FomSystemType, HyperReductionOperatorType>;
  return return_type(trialSpace, fomSystem, hrOp);
}


template<
  class TrialSpaceType,
  class FomSystemType,
  class ResidualMaskerType,
  class JacobianActionMaskerType,
  class HyperReductionOperatorType
  >
auto create_steady_problem(const TrialSpaceType & trialSpace,
			   const FomSystemType & fomSystem,
			   const ResidualMaskerType & rMasker,
			   const JacobianActionMaskerType & jaMasker,
			   const HyperReductionOperatorType & hrOp)
{
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

  // figure out what are the reduced types
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using reduced_residual_type = reduced_state_type;
  using reduced_jac_type = typename
    impl::determine_galerkin_jacobian_type_from_state<reduced_state_type>::type;

  // the hr operator must meet the concept
  static_assert(SteadyGalerkinHyperReducer<
		HyperReductionOperatorType, reduced_residual_type, reduced_jac_type>::value,
		"HyperReductionOperatorType does not meet the SteadyGalerkinHyperReducer");
  // the hr operator must operate on types that are consistent with
  // those deducible from the masker
  static_assert(std::is_same<
		typename HyperReductionOperatorType::residual_operand_type,
		typename ResidualMaskerType::result_type>::value == true,
		"mismatching types");
  static_assert(std::is_same<
		typename HyperReductionOperatorType::jacobian_action_operand_type,
		typename JacobianActionMaskerType::result_type>::value == true,
		"mismatching types");

  using return_type = impl::GalerkinSteadyMaskedSystem<
    reduced_state_type, reduced_residual_type, reduced_jac_type,
    TrialSpaceType, FomSystemType, ResidualMaskerType,
    JacobianActionMaskerType, HyperReductionOperatorType>;
  return return_type(trialSpace, fomSystem, rMasker, jaMasker, hrOp);
}

}}} // end pressio::rom::galerkin
#endif
