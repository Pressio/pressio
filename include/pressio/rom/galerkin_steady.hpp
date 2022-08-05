
#ifndef PRESSIO_ROM_GALERKIN_STEADY_HPP_
#define PRESSIO_ROM_GALERKIN_STEADY_HPP_

#include "impl/reduced_operators_helpers.hpp"
#include "impl/galerkin_steady_default_system.hpp"
#include "impl/galerkin_steady_masked_system.hpp"
#include "impl/galerkin_steady_hypred_system.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

template<
  class TrialSpaceType,
  class FomSystemType,
  mpl::enable_if_t<
    // check for trial concept since affine space subsumes trial concept
    TrialSubspace<TrialSpaceType>::value
    && SteadyFomWithJacobianAction<
      FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
  >
auto create_default_problem(const TrialSpaceType & trialSpace,
			    const FomSystemType & fomObject)
{
  // reduced state and residual have same type
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using reduced_residual_type = reduced_state_type;

  // figure out what is the reduced jacobian type from the state
  using reduced_jac_type = typename impl::determine_galerkin_jacobian_type_from_state<reduced_state_type>::type;

  using return_type = impl::GalerkinSteadyDefaultSystem<
    reduced_state_type, reduced_residual_type,
    reduced_jac_type, TrialSpaceType, FomSystemType>;
  return return_type(trialSpace, fomObject);
}

template<
  class TrialSpaceType,
  class FomSystemType,
  class HyperreductionOperator,
  mpl::enable_if_t<
    // check for trial concept since affine space subsumes trial concept
    TrialSubspace<TrialSpaceType>::value
    && SteadyFomWithJacobianAction<
      FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
  >
auto create_hyperreduced_problem(const TrialSpaceType & trialSpace,
				 const FomSystemType & fomObject,
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
  return return_type(trialSpace, fomObject, hrOp);
}

template<
  class TrialSpaceType,
  class FomSystemType,
  class ResidualMaskerType,
  class JacobianActionMaskerType,
  class HyperreductionOperator,
  mpl::enable_if_t<
    // check for trial concept since affine space subsumes trial concept
    TrialSubspace<TrialSpaceType>::value
    && SteadyFomWithJacobianAction<
      FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
  >
auto create_masked_problem(const TrialSpaceType & trialSpace,
			   const FomSystemType & fomObject,
			   const ResidualMaskerType & rMasker,
			   const JacobianActionMaskerType & jaMasker,
			   const HyperreductionOperator & hrOp)
{

  // reduced state and residual have same type
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using reduced_residual_type = reduced_state_type;

  // figure out what is the reduced jacobian type from the state
  using reduced_jac_type = typename impl::determine_galerkin_jacobian_type_from_state<reduced_state_type>::type;

  using return_type = impl::GalerkinSteadyMaskedSystem<
    reduced_state_type, reduced_residual_type, reduced_jac_type,
    TrialSpaceType, FomSystemType, ResidualMaskerType,
    JacobianActionMaskerType, HyperreductionOperator>;
  return return_type(trialSpace, fomObject, rMasker, jaMasker, hrOp);
}

}}} // end pressio::rom::galerkin
#endif
