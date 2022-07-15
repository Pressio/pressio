
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
auto create_default_problem(TrialSpaceType & trialSpace,
         const FomSystemType & fomObj)
{
  return impl::lspg_steady_create_default_problem(trialSpace, fomObj);
}

template<
  class TrialSpaceType, class FomSystemType, class PreconditionerType,
  mpl::enable_if_t<
    ImplicitTimeInvariantFomWithJacobianAction<
      FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
  >
auto create_preconditioned_problem(TrialSpaceType & trialSpace,
          const FomSystemType & fomObj,
          const PreconditionerType & preconditioner)
{
  return impl::lspg_steady_create_preconditioned_problem(trialSpace, fomObj, preconditioner);
}


}}} // end pressio::rom::galerkin
#endif
