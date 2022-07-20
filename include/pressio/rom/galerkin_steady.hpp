
#ifndef PRESSIO_ROM_GALERKIN_STEADY_HPP_
#define PRESSIO_ROM_GALERKIN_STEADY_HPP_

#include "impl/galerkin_steady_create.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

template<
  class TrialSpaceType,
  class FomSystemType,
  mpl::enable_if_t<
    // check for trial concept since affine space subsumes trial concept
    TrialSubspace<TrialSpaceType>::value
    && ImplicitTimeInvariantFomWithJacobianAction<
      FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
  >
auto create_default_problem(const TrialSpaceType & trialSpaceObject,
			    const FomSystemType & fomObject)
{
  return impl::galerkin_create_default_problem(trialSpaceObject, fomObject);
}

template<
  class TrialSpaceType,
  class FomSystemType,
  class HyperreductionOperator,
  mpl::enable_if_t<
    // check for trial concept since affine space subsumes trial concept
    TrialSubspace<TrialSpaceType>::value
    && ImplicitTimeInvariantFomWithJacobianAction<
      FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
  >
auto create_hyperreduced_problem(const TrialSpaceType & trialSpaceObject,
				 const FomSystemType & fomObject,
				 const HyperreductionOperator & hrOp)
{
  return impl::galerkin_create_hyperreduced_problem(trialSpaceObject, fomObject, hrOp);
}

template<
  class TrialSpaceType,
  class FomSystemType,
  class MaskerType,
  class HyperreductionOperator,
  mpl::enable_if_t<
    // check for trial concept since affine space subsumes trial concept
    TrialSubspace<TrialSpaceType>::value
    && ImplicitTimeInvariantFomWithJacobianAction<
      FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
  >
auto create_masked_problem(const TrialSpaceType & trialSpaceObject,
			   const FomSystemType & fomObject,
			   const MaskerType & masker,
			   const HyperreductionOperator & hrOp)
{
  return impl::galerkin_create_masked_problem(trialSpaceObject, fomObject, masker, hrOp);
}

}}} // end pressio::rom::galerkin
#endif
