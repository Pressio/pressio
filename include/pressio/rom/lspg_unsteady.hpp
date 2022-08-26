
#ifndef PRESSIO_ROM_LSPG_UNSTEADY_HPP_
#define PRESSIO_ROM_LSPG_UNSTEADY_HPP_

#include "./impl/lspg_unsteady_fom_states_manager.hpp"
#include "./impl/lspg_unsteady_default_rj_policy.hpp"
#include "./impl/lspg_unsteady_problem.hpp"

namespace pressio{ namespace rom{ namespace lspg{

template<
  class TrialSpaceType,
  class FomSystemType
>
auto create_default_problem(::pressio::ode::StepScheme schemeName,
			    const TrialSpaceType & trialSpace,
			    const FomSystemType & fomSystem)
{

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

}}} // end pressio::rom::galerkin
#endif
