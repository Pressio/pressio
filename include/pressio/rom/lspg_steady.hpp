
#ifndef PRESSIO_ROM_LSPG_STEADY_HPP_
#define PRESSIO_ROM_LSPG_STEADY_HPP_

#include "impl/lspg_steady_create.hpp"

namespace pressio{ namespace rom{ namespace lspg{

template<
  class TrialSpaceType, class FomSystemType,
  mpl::enable_if_t<
    ImplicitTimeInvariantFomWithJacobianAction<
      FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
  >
auto create_default_problem(TrialSpaceType & trialSpaceObject,
			    const FomSystemType & fomObject)
{
  return impl::lspg_steady_create_default_problem(trialSpaceObject, fomObject);
}

template<
  class TrialSpaceType, class FomSystemType, class PreconditionerType,
  mpl::enable_if_t<
    ImplicitTimeInvariantFomWithJacobianAction<
      FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
  >
auto create_preconditioned_problem(TrialSpaceType & trialSpaceObject,
				   const FomSystemType & fomObject,
				   const PreconditionerType & preconditioner)
{
  return impl::lspg_steady_create_preconditioned_problem(trialSpaceObject, fomObject, preconditioner);
}

template<
  class TrialSpaceType, class FomSystemType,
  mpl::enable_if_t<
    ImplicitTimeInvariantFomWithJacobianAction<
      FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
  >
auto create_hyperreduced_problem(TrialSpaceType & trialSpaceObject,
				 const FomSystemType & fomObject)
{
  // impl-wise, a steady hypred problem is equivalent to a default problem
  return impl::lspg_steady_create_default_problem(trialSpaceObject, fomObject);
}

template<
  class TrialSpaceType, class FomSystemType, class MaskerType,
  mpl::enable_if_t<
    ImplicitTimeInvariantFomWithJacobianAction<
      FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
  >
auto create_masked_problem(TrialSpaceType & trialSpaceObject,
			   const FomSystemType & fomObject,
			   const MaskerType & masker)
{
  return impl::lspg_steady_create_masked_problem(trialSpaceObject,
						 fomObject, masker);
}

}}} // end pressio::rom::galerkin
#endif
