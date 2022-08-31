
#ifndef PRESSIO_ROM_LSPG_UNSTEADY_HPP_
#define PRESSIO_ROM_LSPG_UNSTEADY_HPP_

#include "./impl/lspg_unsteady_fom_states_manager.hpp"
#include "./impl/lspg_unsteady_fully_discrete_system.hpp"
#include "./impl/lspg_unsteady_rj_policy_default.hpp"
#include "./impl/lspg_unsteady_rj_policy_hypred.hpp"
#include "./impl/lspg_unsteady_mask_decorator.hpp"
#include "./impl/lspg_unsteady_problem.hpp"

namespace pressio{ namespace rom{

namespace impl{
void valid_scheme_for_lspg_else_throw(::pressio::ode::StepScheme name){
  if (   name != ::pressio::ode::StepScheme::BDF1
      && name != ::pressio::ode::StepScheme::BDF2)
  {
    throw std::runtime_error("LSPG currently accepting BDF1 or BDF2");
  }
}
}//end impl

namespace lspg{

template<
  std::size_t num_states,
  class TrialSpaceType,
  class FomSystemType
>
auto create_unsteady_problem(const TrialSpaceType & trialSpace,
			     const FomSystemType & fomSystem)
{

  // sufficient to satisfy the TrialSubspace concept since
  // the AffineSpace concept subsumes the TrialSubspace one
  static_assert(TrialSubspace<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialSubspace concept");

  static_assert(FullyDiscreteSystemWithJacobianAction<
		FomSystemType, num_states, typename TrialSpaceType::basis_type>::value,
		"FomSystemType does not meet the FullyDiscreteSystemWithJacobianAction concept");

  static_assert(std::is_same<typename TrialSpaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using lspg_residual_type = typename FomSystemType::discrete_residual_type;
  using lspg_jacobian_type = typename TrialSpaceType::basis_type;

  using system_type = impl::LspgFullyDiscreteSystem<
    num_states, independent_variable_type, reduced_state_type,
    lspg_residual_type, lspg_jacobian_type,
    TrialSpaceType, FomSystemType>;

  using return_type = impl::LspgUnsteadyProblem<TrialSpaceType, system_type, num_states>;
  return return_type(trialSpace, fomSystem);
}






template<
  class TrialSpaceType,
  class FomSystemType
>
auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,
			     const TrialSpaceType & trialSpace,
			     const FomSystemType & fomSystem)
{

  impl::valid_scheme_for_lspg_else_throw(schemeName);

  // sufficient to satisfy the TrialSubspace concept since
  // the AffineSpace concept subsumes the TrialSubspace one
  static_assert(TrialSubspace<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialSubspace concept");

  static_assert(SemiDiscreteFomWithJacobianAction<
		FomSystemType, typename TrialSpaceType::basis_type>::value,
		"FomSystemType does not meet the SemiDiscreteFomWithJacobianAction concept");

  static_assert(std::is_same<typename TrialSpaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;

  /* use the FOM rhs to represent the lspg residual
     and the fom jacobian action on basis as the lspg jacobian */
  using lspg_residual_type = typename FomSystemType::right_hand_side_type;
  using lspg_jacobian_type = typename TrialSpaceType::basis_type;

  // defining an unsteady lspg problem boils down to
  // defining a "custom" residual and jacobian policy
  // since for lspg we need to customize how we do time stepping
  using rj_policy_type = impl::LspgUnsteadyResidualJacobianPolicy<
    independent_variable_type, reduced_state_type,
    lspg_residual_type, lspg_jacobian_type,
    TrialSpaceType, FomSystemType>;

  using return_type = impl::LspgUnsteadyProblem<TrialSpaceType, rj_policy_type>;
  return return_type(schemeName, trialSpace, fomSystem);
}

template<
  class TrialSpaceType,
  class FomSystemType,
  class HypRedUpdaterType
>
auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,
			     const TrialSpaceType & trialSpace,
			     const FomSystemType & fomSystem,
			     const HypRedUpdaterType & hypRedUpdater)
{

  impl::valid_scheme_for_lspg_else_throw(schemeName);

  // sufficient to satisfy the TrialSubspace concept since
  // the AffineSpace concept subsumes the TrialSubspace one
  static_assert(TrialSubspace<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialSubspace concept");

  static_assert(SemiDiscreteFomWithJacobianAction<
		FomSystemType, typename TrialSpaceType::basis_type>::value,
		"FomSystemType does not meet the SemiDiscreteFomWithJacobianAction concept");

  static_assert(std::is_same<typename TrialSpaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  /* use the FOM rhs to represent the lspg residual
     and the fom jacobian action on basis as the lspg jacobian */
  using lspg_residual_type = typename FomSystemType::right_hand_side_type;
  using lspg_jacobian_type = typename TrialSpaceType::basis_type;

  // defining an unsteady lspg problem boils down to
  // defining a "custom" residual and jacobian policy
  // since for lspg we need to customize how we do time stepping
  using rj_policy_type = impl::LspgUnsteadyResidualJacobianPolicyHypRed<
    independent_variable_type, reduced_state_type,
    lspg_residual_type, lspg_jacobian_type,
    TrialSpaceType, FomSystemType, HypRedUpdaterType>;

  using return_type = impl::LspgUnsteadyProblem<TrialSpaceType, rj_policy_type>;
  return return_type(schemeName, trialSpace, fomSystem, hypRedUpdater);
}

template<
  class TrialSpaceType,
  class FomSystemType,
  class ResidualMaskerType,
  class JacobianActionMaskerType
>
auto create_unsteady_problem(::pressio::ode::StepScheme schemeName,
			     const TrialSpaceType & trialSpace,
			     const FomSystemType & fomSystem,
			     const ResidualMaskerType & rMasker,
			     const JacobianActionMaskerType & jaMasker)
{

  impl::valid_scheme_for_lspg_else_throw(schemeName);

  // sufficient to satisfy the TrialSubspace concept since
  // the AffineSpace concept subsumes the TrialSubspace one
  static_assert(TrialSubspace<TrialSpaceType>::value,
		"TrialSpaceType does not meet the TrialSubspace concept");

  static_assert(SemiDiscreteFomWithJacobianAction<
		FomSystemType, typename TrialSpaceType::basis_type>::value,
		"FomSystemType does not meet the SemiDiscreteFomWithJacobianAction concept");

  static_assert(std::is_same<typename TrialSpaceType::full_state_type,
		typename FomSystemType::state_type>::value == true,
		"Mismatching fom states");

  // masker must be invariant in time
  static_assert(TimeInvariantMasker<ResidualMaskerType>::value,
		"ResidualMaskerType does not meet the TimeInvariantMasker concept");
  static_assert(TimeInvariantMasker<JacobianActionMaskerType>::value,
		"JacobianActionMaskerType does not meet the TimeInvariantMasker concept");

  // ensure the masker acts on the FOM types
  static_assert(std::is_same<
		typename ResidualMaskerType::operand_type,
		typename FomSystemType::right_hand_side_type>::value == true,
		"mismatching types of rhs masker and fom right_hand_side_type");
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<typename TrialSpaceType::basis_type const &>()) );
  static_assert(std::is_same<
		typename JacobianActionMaskerType::operand_type,
		fom_jac_action_result_type>::value == true,
		"mismatching types of jacobian action masker and fom");

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;

  using lspg_unmasked_residual_type = typename FomSystemType::right_hand_side_type;
  using lspg_unmasked_jacobian_type = typename TrialSpaceType::basis_type;

  using lspg_residual_type = typename ResidualMaskerType::result_type;
  using lspg_jacobian_type = typename JacobianActionMaskerType::result_type;

  // defining an unsteady lspg problem boils down to
  // defining a "custom" residual and jacobian policy
  // since for lspg we need to customize how we do time stepping
  using rj_policy_type =
    impl::LspgMaskDecorator<
      ResidualMaskerType, JacobianActionMaskerType,
    impl::LspgUnsteadyResidualJacobianPolicy<
      independent_variable_type, reduced_state_type,
      lspg_unmasked_residual_type, lspg_unmasked_jacobian_type,
      TrialSpaceType, FomSystemType
      >
    >;

  using return_type = impl::LspgUnsteadyProblem<TrialSpaceType, rj_policy_type>;
  return return_type(schemeName, trialSpace, fomSystem, rMasker, jaMasker);
}

}}} // end pressio::rom::lspg
#endif
