
#ifndef PRESSIO_ROM_LSPG_STEADY_HPP_
#define PRESSIO_ROM_LSPG_STEADY_HPP_

#include "./impl/reduced_operators_helpers.hpp"
#include "./impl/lspg_steady_default_system.hpp"

namespace pressio{ namespace rom{ namespace lspg{

template<
  class TrialSpaceType,
  class FomSystemType,
  mpl::enable_if_t<
       // sufficient to satisfy the TrialSubspace concept since
       // the AffineSpace concept subsumes the TrialSubspace one
       TrialSubspace<TrialSpaceType>::value
    && SteadyFomWithJacobianAction<
      FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
  >
auto create_default_problem(const TrialSpaceType & trialSpace,
			    const FomSystemType & fomSystem)
{

  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using system_type = impl::LspgSteadyDefaultSystem<
    reduced_state_type, TrialSpaceType, FomSystemType>;
  return system_type(trialSpace, fomSystem);
}

// template<
//   class TrialSpaceType, class FomSystemType,
//   mpl::enable_if_t<
//     SteadyFomWithJacobianAction<
//       FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
//   >
// auto create_hyperreduced_problem(TrialSpaceType & trialSpaceObject,
// 				 const FomSystemType & fomSystem)
// {
//   // impl-wise, it is equivalent to a default problem
//   return impl::lspg_steady_create_default_problem(trialSpaceObject, fomSystem);
// }

// template<
//   class TrialSpaceType, class FomSystemType, class MaskerType,
//   mpl::enable_if_t<
//     SteadyFomWithJacobianAction<
//       FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
//   >
// auto create_masked_problem(TrialSpaceType & trialSpaceObject,
// 			   const FomSystemType & fomSystem,
// 			   const MaskerType & masker)
// {
// using reduced_state_type = typename TrialSpaceType::reduced_state_type;
// using system_type = impl::LspgSteadyMaskedSystem<
//   reduced_state_type, TrialSpaceType, FomSystemType, MaskerType>;
// return system_type(trialSpace, fomObj, masker);
// }

// //
// // overload set for preconditioned problems
// // (might not need these if the preconditioner is instead provided to solver)
// //
// template<
//   class TrialSpaceType, class FomSystemType, class PreconditionerType,
//   mpl::enable_if_t<
//     SteadyFomWithJacobianAction<
//       FomSystemType, typename TrialSpaceType::basis_type>::value, int > = 0
//   >
// auto create_preconditioned_problem(TrialSpaceType & trialSpaceObject,
//            const FomSystemType & fomSystem,
//            const PreconditionerType & preconditioner)
// {
// using reduced_state_type = typename TrialSpaceType::reduced_state_type;
// using system_type = impl::LspgSteadyPreconditionedSystem<
//   reduced_state_type, TrialSpaceType, FomSystemType, PreconditionerType>;
// return system_type(trialSpace, fomObj, preconditioner);
// }

}}} // end pressio::rom::galerkin
#endif
