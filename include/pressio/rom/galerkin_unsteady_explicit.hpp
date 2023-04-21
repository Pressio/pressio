
#ifndef ROM_GALERKIN_UNSTEADY_EXPLICIT_HPP_
#define ROM_GALERKIN_UNSTEADY_EXPLICIT_HPP_

#include "./reduced_operators_traits.hpp"
#include "impl/galerkin_helpers.hpp"
#include "impl/galerkin_unsteady_explicit_problem.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_only.hpp"
#include "impl/galerkin_unsteady_system_hypred_rhs_only.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_with_mass_matrix.hpp"
#include "impl/galerkin_unsteady_system_masked_rhs_only.hpp"

namespace pressio{ namespace rom{ namespace galerkin{


// -------------------------------------------------------------
// default
// -------------------------------------------------------------

#ifdef PRESSIO_ENABLE_CXX20
template<class TrialSubspaceType, class FomSystemType>
  requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
  && RealValuedSemiDiscreteFom<FomSystemType>
  && std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
#else
template<
 class TrialSubspaceType, class FomSystemType,
 mpl::enable_if_t<
  PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value
   && RealValuedSemiDiscreteFom<FomSystemType>::value
   && !RealValuedSemiDiscreteFomWithMassMatrixAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value
   && std::is_same_v<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
   , int> = 0
 >
#endif
auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,  /*(1)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem)
{

  impl::valid_scheme_for_explicit_galerkin_else_throw(schemeName, "galerkin_default_explicit");
  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using reduced_rhs_type = impl::explicit_galerkin_default_reduced_rhs_t<TrialSubspaceType>;
  // the "system" implements the math
  using galerkin_system = impl::GalerkinDefaultOdeSystemOnlyRhs<
    ind_var_type, reduced_state_type, reduced_rhs_type, TrialSubspaceType, FomSystemType>;

  // a Galerkin problem contains (beside other things) a pressio stepper.
  // the reason for this is that a problem can potentially expose
  // more methods than what the underlying stepper does.
  using return_type = impl::GalerkinUnsteadyExplicitProblem<galerkin_system>;
  return return_type(schemeName, trialSpace, fomSystem);
}

// -------------------------------------------------------------
// default with mass matrix
// -------------------------------------------------------------
#ifdef PRESSIO_ENABLE_CXX20
template<class TrialSubspaceType, class FomSystemType>
  requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
  && RealValuedSemiDiscreteFomWithMassMatrixAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
  && std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
#else
template<
  class TrialSubspaceType, class FomSystemType,
  mpl::enable_if_t<
    PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value
    && RealValuedSemiDiscreteFomWithMassMatrixAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value
    && std::is_same_v<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
    , int > = 0
  >
#endif
auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,  /*(2)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem)
{

  impl::valid_scheme_for_explicit_galerkin_else_throw(schemeName, "galerkin_default_explicit");
  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using reduced_rhs_type = impl::explicit_galerkin_default_reduced_rhs_t<TrialSubspaceType>;
  using reduced_mm_type = impl::explicit_galerkin_default_reduced_mass_matrix_t<TrialSubspaceType>;
  using galerkin_system = impl::GalerkinDefaultOdeSystemOnlyRhsAndMassMatrix<
    ind_var_type, reduced_state_type, reduced_rhs_type,
    reduced_mm_type, TrialSubspaceType, FomSystemType>;

  using return_type = impl::GalerkinUnsteadyWithMassMatrixExplicitProblem<galerkin_system>;
  return return_type(schemeName, trialSpace, fomSystem);
}

// -------------------------------------------------------------
// hyper-reduced
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType>
#ifdef PRESSIO_ENABLE_CXX20
requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
&& RealValuedSemiDiscreteFom<FomSystemType>
&& std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
#endif
auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,  /*(3)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const HyperReducerType & hyperReducer)
{

  impl::valid_scheme_for_explicit_galerkin_else_throw(schemeName, "galerkin_hypred_explicit");

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using reduced_rhs_type = impl::explicit_galerkin_default_reduced_rhs_t<TrialSubspaceType>;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinHyperReducedOdeSystemOnlyRhs<
    ind_var_type, reduced_state_type, reduced_rhs_type,
    TrialSubspaceType, FomSystemType, HyperReducerType>;

  using return_type = impl::GalerkinUnsteadyExplicitProblem<galerkin_system>;
  return return_type(schemeName, trialSpace, fomSystem, hyperReducer);
}

// -------------------------------------------------------------
// masked
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class HyperReducerType>
#ifdef PRESSIO_ENABLE_CXX20
requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
&& RealValuedSemiDiscreteFom<FomSystemType>
&& std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
#endif
auto create_unsteady_explicit_problem(::pressio::ode::StepScheme schemeName,  /*(4)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const MaskerType & masker,
				      const HyperReducerType & hyperReducer)
{

  impl::valid_scheme_for_explicit_galerkin_else_throw(schemeName, "galerkin_masked_explicit");
  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using reduced_rhs_type = reduced_state_type;
  // the "system" implements the math
  using galerkin_system = impl::GalerkinMaskedOdeSystemOnlyRhs<
    ind_var_type, reduced_state_type, reduced_rhs_type,
    TrialSubspaceType, FomSystemType, MaskerType, HyperReducerType>;

  using return_type = impl::GalerkinUnsteadyExplicitProblem<galerkin_system>;
  return return_type(schemeName, trialSpace, fomSystem, masker, hyperReducer);
}

}}} // end pressio::rom::galerkin
#endif  // ROM_GALERKIN_UNSTEADY_EXPLICIT_HPP_
