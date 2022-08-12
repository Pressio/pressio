
#ifndef PRESSIO_ROM_GALERKIN_UNSTEADY_HPP_
#define PRESSIO_ROM_GALERKIN_UNSTEADY_HPP_

#include "impl/reduced_operators_helpers.hpp"
#include "impl/galerkin_unsteady_explicit_problem.hpp"
#include "impl/galerkin_unsteady_systems.hpp"

namespace pressio{ namespace rom{

namespace impl{
void explicit_scheme_else_throw(::pressio::ode::StepScheme name,
				const std::string & str)
{
  if (!::pressio::ode::is_explicit_scheme(name)){
    throw std::runtime_error(str + " requires an explicit stepper");
  }
}

void implicit_scheme_else_throw(::pressio::ode::StepScheme name,
				const std::string & str)
{
  if (!::pressio::ode::is_implicit_scheme(name)){
    throw std::runtime_error(str + " requires an implicit stepper");
  }
}
}//end impl

namespace galerkin{

template<
  class TrialSpaceType,
  class FomSystemType,
  mpl::enable_if_t<
       // sufficient to satisfy the TrialSubspace concept since
       // the AffineSpace concept subsumes the TrialSubspace one
       TrialSubspace<TrialSpaceType>::value
       // this overload is for the case WITHOUT mass matrix
    && SemiDiscreteFom<FomSystemType>::value
    && !SemiDiscreteFomWithMassMatrixAction<
      FomSystemType, typename TrialSpaceType::basis_type >::value
    && !SemiDiscreteFomWithConstantMassMatrixAction<
      FomSystemType, typename TrialSpaceType::basis_type >::value,
    int > = 0
  >
auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,
				      const TrialSpaceType & trialSpace,
				      const FomSystemType & fomSystem)
{
  impl::explicit_scheme_else_throw(schemeName, "galerkin_default_explicit");

  using independent_variable_type = typename FomSystemType::time_type;

  // for the reduced rhs, use the type of the reduced state
  using reduced_state_type = typename TrialSpaceType::reduced_state_type;
  using reduced_rhs_type = reduced_state_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinDefaultOdeSystemOnlyRhs<
    independent_variable_type, reduced_state_type, reduced_rhs_type,
    TrialSpaceType, FomSystemType>;

  // a Galerkin problem is simply a wrapper of a pressio stepper.
  // the reason for this is that a problem can potentially expose
  // more methods than what the underlying stepper does.
  using return_type = impl::GalerkinUnsteadyExplicitProblem<false, galerkin_system>;
  return return_type(schemeName, trialSpace, fomSystem);
}

}}} // end pressio::rom::galerkin
#endif

// //
// // not yet finalizd
// //
// template<
//   class TrialSpaceType,
//   class FomSystemType,
//   mpl::enable_if_t<
//     // check for trial concept since affine space subsumes trial concept
//     TrialSubspace<TrialSpaceType>::value
//     && SemiDiscreteFomWithRhsAndMassMatrixAction<
//       FomSystemType, typename TrialSpaceType::basis_type >::value, int > = 0
//   >
// auto create_default_explicit_problem(::pressio::ode::StepScheme schemeName,
// 				     const TrialSpaceType & trialSpaceObject,
// 				     const FomSystemType & fomSystem)
// {

//   impl::explicit_scheme_else_throw(schemeName, "galerkin_default_explicit");
//   using independent_variable_type = typename FomSystemType::time_type;

//   // rom state and rhs have same type
//   using reduced_state_type = typename TrialSpaceType::reduced_state_type;
//   using reduced_rhs_type = reduced_state_type;

//   // figure out what is the mass matrix type from the state
//   using rom_mm_type = typename impl::galerkin_mass_matrix_type_if_state_is<
//     reduced_state_type>::type;

//   using galerkin_system = impl::GalerkinUnsteadyDefaultOdeSystemWithMassMatrix<
//     independent_variable_type, reduced_state_type, reduced_rhs_type, rom_mm_type,
//     TrialSpaceType, FomSystemType>;

//   using return_type = impl::GalerkinUnsteadyExplicitProblem<true, galerkin_system>;
//   return return_type(schemeName, trialSpace, fomObj);
// }

// template<
//   class TrialSpaceType,
//   class FomSystemType,
//   class HyperreductionOperator,
//   mpl::enable_if_t<
//     // check for trial concept since affine space subsumes trial concept
//     TrialSubspace<TrialSpaceType>::value
//     // constraint for case WITHOUT mass matrix
//     && SemiDiscreteFomWithRhs<FomSystemType>::value
//     && !SemiDiscreteFomWithRhsAndMassMatrixAction<
//       FomSystemType, typename TrialSpaceType::basis_type >::value, int > = 0
//   >
// auto create_hyperreduced_explicit_problem(::pressio::ode::StepScheme schemeName,
// 					  const TrialSpaceType & trialSpaceObject,
// 					  const FomSystemType & fomSystem,
// 					  const HyperreductionOperator & hrOp)
// {

//   impl::explicit_scheme_else_throw(schemeName, "galerkin_default_explicit");

//   using independent_variable_type = typename FomSystemType::time_type;

//   // rom state and rhs have same type
//   using reduced_state_type = typename TrialSpaceType::reduced_state_type;
//   using reduced_rhs_type = reduced_state_type;

//   using galerkin_system = impl::GalerkinUnsteadyHypRedOdeSystem<
//     independent_variable_type, reduced_state_type, reduced_rhs_type,
//     TrialSpaceType, FomSystemType, HyperreductionOperator>;

//   using return_type = impl::GalerkinUnsteadyExplicitProblem<false, galerkin_system>;
//   return return_type(schemeName, trialSpace, fomObj, hrOp);
// }

// template<
//   class TrialSpaceType,
//   class FomSystemType,
//   class HyperreductionOperator,
//   mpl::enable_if_t<
//     // check for trial concept since affine space subsumes trial concept
//     TrialSubspace<TrialSpaceType>::value
//     && SemiDiscreteFomWithRhsAndConstantMassMatrixAction<
//       FomSystemType, typename TrialSpaceType::basis_type >::value, int > = 0
//   >
// auto create_hyperreduced_explicit_problem(::pressio::ode::StepScheme schemeName,
// 					  const TrialSpaceType & trialSpaceObject,
// 					  const FomSystemType & fomSystem,
// 					  const HyperreductionOperator & hrOp)
// {

//   impl::explicit_scheme_else_throw(schemeName, "galerkin_default_explicit");

//   using independent_variable_type = typename FomSystemType::time_type;

//   // rom state and rhs have same type
//   using reduced_state_type = typename TrialSpaceType::reduced_state_type;
//   using reduced_rhs_type = reduced_state_type;

//   // figure out what is the mass matrix type from the state
//   using rom_mm_type = typename impl::galerkin_mass_matrix_type_if_state_is<
//     reduced_state_type>::type;

//   using galerkin_system = impl::GalerkinUnsteadyHypRedOdeSystemWithConstantMassMatrix<
//     independent_variable_type, reduced_state_type, reduced_rhs_type, rom_mm_type,
//     TrialSpaceType, FomSystemType, HyperreductionOperator>;

//   using return_type = impl::GalerkinUnsteadyExplicitProblem<true, galerkin_system>;
//   return return_type(schemeName, trialSpace, fomObj, hrOp);
// }

// template<
//   class TrialSpaceType,
//   class FomSystemType,
//   class MaskerType,
//   class HyperreductionOperator,
//   mpl::enable_if_t<
//     // check for trial concept since affine space subsumes trial concept
//     TrialSubspace<TrialSpaceType>::value
//     // constraint for case WITHOUT mass matrix
//     && SemiDiscreteFomWithRhs<FomSystemType>::value
//     && !SemiDiscreteFomWithRhsAndMassMatrixAction<
//       FomSystemType, typename TrialSpaceType::basis_type >::value, int > = 0
//   >
// auto create_masked_explicit_problem(::pressio::ode::StepScheme schemeName,
// 				    const TrialSpaceType & trialSpaceObject,
// 				    const FomSystemType & fomSystem,
// 				    const MaskerType & masker,
// 				    const HyperreductionOperator & hrOp)
// {

//   impl::explicit_scheme_else_throw(schemeName, "galerkin_default_explicit");

//   using independent_variable_type = typename FomSystemType::time_type;

//   // rom state and rhs have same type
//   using reduced_state_type = typename TrialSpaceType::reduced_state_type;
//   using reduced_rhs_type = reduced_state_type;

//   using galerkin_system = impl::GalerkinUnsteadyMaskedOdeSystem<
//     independent_variable_type, reduced_state_type, reduced_rhs_type,
//     TrialSpaceType, FomSystemType, MaskerType, HyperreductionOperator>;

//   using return_type = impl::GalerkinUnsteadyExplicitProblem<false, galerkin_system>;
//   return return_type(schemeName, trialSpace, fomObj, masker, hrOp);
// }

// template<
//   class TrialSpaceType,
//   class FomSystemType,
//   class MaskerType,
//   class HyperreductionOperator,
//   mpl::enable_if_t<
//     // check for trial concept since affine space subsumes trial concept
//     TrialSubspace<TrialSpaceType>::value
//     && SemiDiscreteFomWithRhsAndConstantMassMatrixAction<
//       FomSystemType, typename TrialSpaceType::basis_type >::value, int > = 0
//   >
// auto create_masked_explicit_problem(::pressio::ode::StepScheme schemeName,
// 				    const TrialSpaceType & trialSpaceObject,
// 				    const FomSystemType & fomSystem,
// 				    const MaskerType & masker,
// 				    const HyperreductionOperator & hrOp)
// {

//   impl::explicit_scheme_else_throw(schemeName, "galerkin_default_explicit");

//   using independent_variable_type = typename FomSystemType::time_type;

//   // rom state and rhs have same type
//   using reduced_state_type = typename TrialSpaceType::reduced_state_type;
//   using reduced_rhs_type = reduced_state_type;

//   // figure out what is the mass matrix type from the state
//   using rom_mm_type = typename impl::galerkin_mass_matrix_type_if_state_is<
//     reduced_state_type>::type;

//   using galerkin_system = impl::GalerkinUnsteadyMaskedOdeSystemWithConstantMassMatrix<
//     independent_variable_type, reduced_state_type, reduced_rhs_type, rom_mm_type,
//     TrialSpaceType, FomSystemType, MaskerType, HyperreductionOperator>;

//   using return_type = impl::GalerkinUnsteadyExplicitProblem<true, galerkin_system>;
//   return return_type(schemeName, trialSpace, fomObj, hrOp);
// }
