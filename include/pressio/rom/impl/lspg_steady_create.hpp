
#ifndef PRESSIO_ROM_LSPG_STEADY_IMPL_CREATE_FUNCTIONS_HPP_
#define PRESSIO_ROM_LSPG_STEADY_IMPL_CREATE_FUNCTIONS_HPP_

#include "reduced_operators_helpers.hpp"
#include "lspg_steady_default_system.hpp"
#include "lspg_steady_masked_system.hpp"
#include "lspg_steady_preconditioned_system.hpp"

namespace pressio{ namespace rom{ namespace impl{

template<
  class TrialSpaceType, class FomSystemType,
  mpl::enable_if_t<
    SteadyFomWithJacobianAction<
      FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
  >
auto lspg_steady_create_default_problem(TrialSpaceType & trialSpace,
					const FomSystemType & fomObj)
{

  using lspg_state_type = typename TrialSpaceType::reduced_state_type;

  using system_type = impl::LspgSteadyDefaultSystem<
    lspg_state_type, TrialSpaceType, FomSystemType>;
  return system_type(trialSpace, fomObj);
}

template<
  class TrialSpaceType, class FomSystemType, class PreconditionerType,
  mpl::enable_if_t<
    SteadyFomWithJacobianAction<
      FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
  >
auto lspg_steady_create_preconditioned_problem(TrialSpaceType & trialSpace,
					       const FomSystemType & fomObj,
					       const PreconditionerType & preconditioner)
{

  using lspg_state_type = typename TrialSpaceType::reduced_state_type;

  using system_type = impl::LspgSteadyPreconditionedSystem<
    lspg_state_type, TrialSpaceType, FomSystemType, PreconditionerType>;
  return system_type(trialSpace, fomObj, preconditioner);
}

template<
  class TrialSpaceType, class FomSystemType, class MaskerType,
  mpl::enable_if_t<
    SteadyFomWithJacobianAction<
      FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
  >
auto lspg_steady_create_masked_problem(TrialSpaceType & trialSpace,
				       const FomSystemType & fomObj,
				       const MaskerType & masker)
{

  using lspg_state_type = typename TrialSpaceType::reduced_state_type;

  using system_type = impl::LspgSteadyMaskedSystem<
    lspg_state_type, TrialSpaceType, FomSystemType, MaskerType>;
  return system_type(trialSpace, fomObj, masker);
}

}}} // end pressio::rom::impl
#endif
