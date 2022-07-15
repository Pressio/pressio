
#ifndef PRESSIO_ROM_GALERKIN_STEADY_IMPL_CREATE_FUNCTIONS_HPP_
#define PRESSIO_ROM_GALERKIN_STEADY_IMPL_CREATE_FUNCTIONS_HPP_

#include "reduced_operators_helpers.hpp"
#include "galerkin_steady_default_system.hpp"
#include "galerkin_steady_masked_system.hpp"
#include "galerkin_steady_hypred_system.hpp"

namespace pressio{ namespace rom{ namespace impl{

template<
  class TrialSpaceType,
  class FomSystemType>
auto galerkin_create_default_problem(const TrialSpaceType & trialSpace,
          const FomSystemType & fomObj)
{
  // reduced state and residual have same type
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using reduced_residual_type = reduced_state_type;

  // figure out what is the reduced jacobian type from the state
  using reduced_jac_type = typename impl::determine_galerkin_jacobian_type_from_state<reduced_state_type>::type;

  using return_type = impl::GalerkinSteadyDefaultSystem<
    reduced_state_type, reduced_residual_type,
    reduced_jac_type, TrialSpaceType, FomSystemType>;
  return return_type(trialSpace, fomObj);
}

template<
  class TrialSpaceType,
  class FomSystemType,
  class HyperreductionOperator>
auto galerkin_create_hyperreduced_problem(const TrialSpaceType & trialSpace,
         const FomSystemType & fomObj,
         const HyperreductionOperator & hrOp)
{
  // reduced state and residual have same type
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using reduced_residual_type = reduced_state_type;

  // figure out what is the reduced jacobian type from the state
  using reduced_jac_type = typename impl::determine_galerkin_jacobian_type_from_state<reduced_state_type>::type;

  using return_type = impl::GalerkinSteadyHypRedSystem<
    reduced_state_type, reduced_residual_type, reduced_jac_type,
    TrialSpaceType, FomSystemType, HyperreductionOperator>;
  return return_type(trialSpace, fomObj, hrOp);
}

template<
  class TrialSpaceType,
  class FomSystemType,
  class MaskerType,
  class HyperreductionOperator
  >
auto galerkin_create_masked_problem(const TrialSpaceType & trialSpace,
         const FomSystemType & fomObj,
         const MaskerType & masker,
         const HyperreductionOperator & hrOp)
{
  // reduced state and residual have same type
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using reduced_residual_type = reduced_state_type;

  // figure out what is the reduced jacobian type from the state
  using reduced_jac_type = typename impl::determine_galerkin_jacobian_type_from_state<reduced_state_type>::type;

  using return_type = impl::GalerkinSteadyMaskedSystem<
    reduced_state_type, reduced_residual_type, reduced_jac_type,
    TrialSpaceType, FomSystemType, MaskerType, HyperreductionOperator>;
  return return_type(trialSpace, fomObj, masker, hrOp);
}


}}} // end pressio::rom::impl
#endif
