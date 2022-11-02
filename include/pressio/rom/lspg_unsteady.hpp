
#ifndef ROM_LSPG_UNSTEADY_HPP_
#define ROM_LSPG_UNSTEADY_HPP_

#include "./impl/lspg_helpers.hpp"
#include "./impl/lspg_unsteady_fom_states_manager.hpp"
#include "./impl/lspg_unsteady_fully_discrete_system.hpp"
#include "./impl/lspg_unsteady_rj_policy_default.hpp"
#include "./impl/lspg_unsteady_rj_policy_hypred.hpp"
#include "./impl/lspg_unsteady_mask_decorator.hpp"
#include "./impl/lspg_unsteady_preconditioning_decorator.hpp"
#include "./impl/lspg_unsteady_problem.hpp"
#if defined PRESSIO_ENABLE_TPL_TRILINOS
#include "./rom_lspg_unsteady_hypred_updater_trilinos.hpp"
#endif

namespace pressio{ namespace rom{ namespace lspg{

// -------------------------------------------------------------
// default
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType>
#ifdef PRESSIO_ENABLE_CXX20
  requires unsteady::ComposableIntoDefaultProblem<TrialSubspaceType, FomSystemType>
#endif
auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,
			     const TrialSubspaceType & trialSpace,
			     const FomSystemType & fomSystem)
{

  impl::valid_scheme_for_lspg_else_throw(schemeName);

#if not defined PRESSIO_ENABLE_CXX20
    static_assert(unsteady::ComposableIntoDefaultProblem<TrialSubspaceType, FomSystemType>::value,
		"does not meet concept");
#endif

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using lspg_residual_type = typename FomSystemType::right_hand_side_type;
  using lspg_jacobian_type =
    decltype(
	     fomSystem.createResultOfJacobianActionOn(trialSpace.basisOfTranslatedSpace())
	     );

  // defining an unsteady lspg problem boils down to
  // defining a "custom" residual and jacobian policy
  // since for lspg we need to customize how we do time stepping
  using rj_policy_type = impl::LspgUnsteadyResidualJacobianPolicy<
    independent_variable_type, reduced_state_type,
    lspg_residual_type, lspg_jacobian_type,
    TrialSubspaceType, FomSystemType>;

  using return_type = impl::LspgUnsteadyProblemSemiDiscreteAPI<TrialSubspaceType, rj_policy_type>;
  return return_type(schemeName, trialSpace, fomSystem);
}

namespace experimental{
// -------------------------------------------------------------
// default with preconditioning
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class PreconditionerType>
#ifdef PRESSIO_ENABLE_CXX20
  requires unsteady::ComposableIntoDefaultProblem<TrialSubspaceType, FomSystemType>
#endif
auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,
			     const TrialSubspaceType & trialSpace,
			     const FomSystemType & fomSystem,
			     const PreconditionerType & preconditioner)
{

  impl::valid_scheme_for_lspg_else_throw(schemeName);

#if not defined PRESSIO_ENABLE_CXX20
    static_assert(unsteady::ComposableIntoDefaultProblem<TrialSubspaceType, FomSystemType>::value,
		"does not meet concept");
#endif

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using lspg_residual_type = typename FomSystemType::right_hand_side_type;
  using lspg_jacobian_type =
    decltype(
	     fomSystem.createResultOfJacobianActionOn(trialSpace.basisOfTranslatedSpace())
	     );

  using rj_policy_type =
    impl::LspgPreconditioningDecorator<
      PreconditionerType, lspg_residual_type, lspg_jacobian_type, TrialSubspaceType,
      impl::LspgUnsteadyResidualJacobianPolicy<
	independent_variable_type, reduced_state_type, lspg_residual_type,
	lspg_jacobian_type, TrialSubspaceType, FomSystemType
	>
    >;

  using return_type = impl::LspgUnsteadyProblemSemiDiscreteAPI<TrialSubspaceType, rj_policy_type>;
  return return_type(schemeName, trialSpace, fomSystem, preconditioner);
}
} //end namespace experimental

// -------------------------------------------------------------
// hyp-red
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class HypRedUpdaterType
#if not defined PRESSIO_ENABLE_CXX20
  , mpl::enable_if_t<
      unsteady::ComposableIntoHyperReducedProblem<TrialSubspaceType, FomSystemType, HypRedUpdaterType>::value,
      int> = 0
#endif
  >
#ifdef PRESSIO_ENABLE_CXX20
  requires unsteady::ComposableIntoHyperReducedProblem<
	TrialSubspaceType, FomSystemType, HypRedUpdaterType>
#endif
auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,
			     const TrialSubspaceType & trialSpace,
			     const FomSystemType & fomSystem,
			     const HypRedUpdaterType & hypRedUpdater)
{

  impl::valid_scheme_for_lspg_else_throw(schemeName);

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using lspg_residual_type = typename FomSystemType::right_hand_side_type;
  using lspg_jacobian_type =
    decltype(
	     fomSystem.createResultOfJacobianActionOn(trialSpace.basisOfTranslatedSpace())
	     );

  using rj_policy_type = impl::LspgUnsteadyResidualJacobianPolicyHypRed<
    independent_variable_type, reduced_state_type,
    lspg_residual_type, lspg_jacobian_type,
    TrialSubspaceType, FomSystemType, HypRedUpdaterType>;

  using return_type = impl::LspgUnsteadyProblemSemiDiscreteAPI<TrialSubspaceType, rj_policy_type>;
  return return_type(schemeName, trialSpace, fomSystem, hypRedUpdater);
}

namespace experimental{
// -------------------------------------------------------------
// hyper-reduced with preconditioning
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class HypRedUpdaterType,
  class PreconditionerType
#if not defined PRESSIO_ENABLE_CXX20
  , mpl::enable_if_t<
      unsteady::ComposableIntoHyperReducedProblem<TrialSubspaceType, FomSystemType, HypRedUpdaterType>::value,
      int> = 0
#endif
  >
#ifdef PRESSIO_ENABLE_CXX20
  requires unsteady::ComposableIntoHyperReducedProblem<
	TrialSubspaceType, FomSystemType, HypRedUpdaterType>
#endif
auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,
			     const TrialSubspaceType & trialSpace,
			     const FomSystemType & fomSystem,
			     const HypRedUpdaterType & hypRedUpdater,
			     const PreconditionerType & preconditioner)
{

  impl::valid_scheme_for_lspg_else_throw(schemeName);

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using lspg_residual_type = typename FomSystemType::right_hand_side_type;
  using lspg_jacobian_type =
    decltype(
	     fomSystem.createResultOfJacobianActionOn(trialSpace.basisOfTranslatedSpace())
	     );

  using rj_policy_type =
    impl::LspgPreconditioningDecorator<
      PreconditionerType, lspg_residual_type, lspg_jacobian_type, TrialSubspaceType,
      impl::LspgUnsteadyResidualJacobianPolicyHypRed<
	independent_variable_type, reduced_state_type, lspg_residual_type,
	lspg_jacobian_type, TrialSubspaceType, FomSystemType, HypRedUpdaterType
	>
    >;

  using return_type = impl::LspgUnsteadyProblemSemiDiscreteAPI<TrialSubspaceType, rj_policy_type>;
  return return_type(schemeName, trialSpace, fomSystem, preconditioner, hypRedUpdater);
}
} //end namespace experimental

// -------------------------------------------------------------
// maskd
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType
#if not defined PRESSIO_ENABLE_CXX20
  , mpl::enable_if_t<
      unsteady::ComposableIntoMaskedProblem<TrialSubspaceType, FomSystemType, MaskerType>::value, int
    > = 0
#endif
  >
#ifdef PRESSIO_ENABLE_CXX20
  requires unsteady::ComposableIntoMaskedProblem<TrialSubspaceType, FomSystemType, MaskerType>
#endif
auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,
			     const TrialSubspaceType & trialSpace,
			     const FomSystemType & fomSystem,
			     const MaskerType & masker)
{

  impl::valid_scheme_for_lspg_else_throw(schemeName);

  using independent_variable_type   = typename FomSystemType::time_type;
  using reduced_state_type          = typename TrialSubspaceType::reduced_state_type;
  using lspg_unmasked_residual_type = typename FomSystemType::right_hand_side_type;
  using lspg_unmasked_jacobian_type = typename TrialSubspaceType::basis_matrix_type;
  using lspg_residual_type =
    decltype(std::declval<MaskerType const>().createResultOfMaskActionOn
	     (std::declval<lspg_unmasked_residual_type const &>()));

  using lspg_jacobian_type =
    decltype(std::declval<MaskerType const>().createResultOfMaskActionOn
	     (std::declval<lspg_unmasked_jacobian_type const &>()));

  using rj_policy_type =
    impl::LspgMaskDecorator<
      MaskerType, lspg_residual_type, lspg_jacobian_type,
      impl::LspgUnsteadyResidualJacobianPolicy<
	independent_variable_type, reduced_state_type,
	lspg_unmasked_residual_type, lspg_unmasked_jacobian_type,
	TrialSubspaceType, FomSystemType
	>
    >;

  using return_type = impl::LspgUnsteadyProblemSemiDiscreteAPI<TrialSubspaceType, rj_policy_type>;
  return return_type(schemeName, trialSpace, fomSystem, masker);
}

// -------------------------------------------------------------
// fully-discrete
// -------------------------------------------------------------
template<
  std::size_t TotalNumberOfDesiredStates,
  class TrialSubspaceType,
  class FomSystemType>
#ifdef PRESSIO_ENABLE_CXX20
  requires FullyDiscreteSystemWithJacobianAction<
	FomSystemType, TotalNumberOfDesiredStates, typename TrialSubspaceType::basis_matrix_type>
#endif
auto create_unsteady_problem(const TrialSubspaceType & trialSpace,
			     const FomSystemType & fomSystem)
{

#if not defined PRESSIO_ENABLE_CXX20
  static_assert(FullyDiscreteSystemWithJacobianAction<
		FomSystemType, TotalNumberOfDesiredStates,
		typename TrialSubspaceType::basis_matrix_type>::value,
		"concept not satisfied");
#endif

  static_assert(std::is_same<typename TrialSubspaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using lspg_residual_type = typename FomSystemType::discrete_residual_type;
  using lspg_jacobian_type =
    decltype(
	     fomSystem.createResultOfDiscreteTimeJacobianActionOn
	     (trialSpace.basisOfTranslatedSpace())
	     );

  using system_type = impl::LspgFullyDiscreteSystem<
    TotalNumberOfDesiredStates, independent_variable_type, reduced_state_type,
    lspg_residual_type, lspg_jacobian_type,
    TrialSubspaceType, FomSystemType>;

  using return_type = impl::LspgUnsteadyProblemFullyDiscreteAPI<
    TotalNumberOfDesiredStates, TrialSubspaceType, system_type>;
  return return_type(trialSpace, fomSystem);
}

}}} // end pressio::rom::lspg
#endif  // ROM_LSPG_UNSTEADY_HPP_
