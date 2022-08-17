
#ifndef PRESSIO_ROM_LSPG_UNSTEADY_HPP_
#define PRESSIO_ROM_LSPG_UNSTEADY_HPP_

#include "./impl/lspg_unsteady_fom_states_manager.hpp"
#include "./impl/lspg_unsteady_default_policy_residual.hpp"
#include "./impl/lspg_unsteady_default_policy_jacobian.hpp"
#include "./impl/lspg_unsteady_problem.hpp"

namespace pressio{ namespace rom{ namespace lspg{

template<
  class TrialSpaceType,
  class FomSystemType,
  mpl::enable_if_t<
       // sufficient to satisfy the TrialSubspace concept since
       // the AffineSpace concept subsumes the TrialSubspace one
       TrialSubspace<TrialSpaceType>::value
       // we need the jacobian action
    && SemiDiscreteFomWithJacobianAction<
	 FomSystemType, typename TrialSpaceType::basis_type>::value
       // this overload is only for jacobian action WITHOUT mass matrix
       // so we need to exclude the concepts with a mass matrix
    && !SemiDiscreteFomComplete<
	 FomSystemType, typename TrialSpaceType::basis_type>::value,
    int > = 0
  >
auto create_default_problem(::pressio::ode::StepScheme schemeName,
			    const TrialSpaceType & trialSpace,
			    const FomSystemType & fomSystem)
{

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;

  /* use the FOM rhs to represent the lspg residual
     and the fom jacobian action on basis as the lspg jacobian
   */
  using lspg_residual_type = typename FomSystemType::right_hand_side_type;
  using lspg_jacobian_type = typename TrialSpaceType::basis_type;

  // defining an unsteady lspg problem boils down to
  // defining a "custom" residual and jacobian policy
  // since for lspg we need to customize how we do time stepping
  using residual_policy_type = impl::LspgUnsteadyResidualPolicy<
    independent_variable_type, reduced_state_type,
    lspg_residual_type, TrialSpaceType, FomSystemType>;
  using jacobian_policy_type = impl::LspgUnsteadyJacobianPolicy<
    independent_variable_type, reduced_state_type,
    lspg_jacobian_type, TrialSpaceType, FomSystemType>;

  using return_type = impl::LspgUnsteadyProblem<
    TrialSpaceType, residual_policy_type, jacobian_policy_type>;
  return return_type(schemeName, trialSpace, fomSystem);
}

}}} // end pressio::rom::galerkin
#endif
