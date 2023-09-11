
#ifndef ROM_GALERKIN_UNSTEADY_IMPLICIT_HPP_
#define ROM_GALERKIN_UNSTEADY_IMPLICIT_HPP_

#include "impl/galerkin_helpers.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_default_rhs_and_jacobian_and_mm.hpp"
#include "impl/galerkin_unsteady_system_hypred_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_masked_rhs_and_jacobian.hpp"
#include "impl/galerkin_unsteady_system_fully_discrete_fom.hpp"
#include "impl/galerkin_unsteady_system_hypred_fully_discrete_fom.hpp"
#include "impl/mixed_fom_rom_unsteady_problem.hpp"
#include "impl/galerkin_unsteady_default_problem_mixed_fom.hpp"

namespace pressio{ namespace rom{ namespace galerkin{

// -------------------------------------------------------------
// default
// -------------------------------------------------------------

#ifdef PRESSIO_ENABLE_CXX20
template<class TrialSubspaceType, class FomSystemType>
  requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
  && RealValuedSemiDiscreteFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
  && std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
#else
template<
 class TrialSubspaceType, class FomSystemType,
 mpl::enable_if_t<
   PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value
   && RealValuedSemiDiscreteFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value
   && !RealValuedSemiDiscreteFomWithJacobianAndMassMatrixAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value
   && std::is_same<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>::value
   , int > = 0
  >
#endif
auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,   /*(1)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem)
{

  impl::valid_scheme_for_implicit_galerkin_else_throw(schemeName, "galerkin_default_implicit");

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types      = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinDefaultOdeSystemRhsAndJacobian<
    ind_var_type, reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType>;

  galerkin_system galSystem(trialSpace, fomSystem);
  return ::pressio::ode::create_implicit_stepper(schemeName, std::move(galSystem));
}

// -------------------------------------------------------------
// default with mass matrix
// -------------------------------------------------------------

#ifdef PRESSIO_ENABLE_CXX20
template<class TrialSubspaceType, class FomSystemType>
  requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
  && RealValuedSemiDiscreteFomWithJacobianAndMassMatrixAction<
      FomSystemType, typename TrialSubspaceType::basis_matrix_type>
  && std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
#else
template<
  class TrialSubspaceType, class FomSystemType,
  mpl::enable_if_t<
    PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>::value
    && RealValuedSemiDiscreteFomWithJacobianAndMassMatrixAction<
         FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value
    && std::is_same<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>::value
    , int > = 0
  >
#endif
auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,   /*(2)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem)
{

  impl::valid_scheme_for_implicit_galerkin_else_throw(schemeName, "galerkin_default_implicit");

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types      = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;
  using reduced_mm_type       = typename default_types::reduced_mass_matrix_type;

  using galerkin_system = impl::GalerkinDefaultOdeSystemRhsJacobianMassMatrix<
    ind_var_type, reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, reduced_mm_type, TrialSubspaceType, FomSystemType>;

  galerkin_system galSystem(trialSpace, fomSystem);
  return ::pressio::ode::create_implicit_stepper(schemeName, std::move(galSystem));
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
&& RealValuedSemiDiscreteFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
&& std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
#endif
auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,   /*(3)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const HyperReducerType & hyperReducer)
{

  impl::valid_scheme_for_implicit_galerkin_else_throw(schemeName, "galerkin_hypred_implicit");

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types      = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinHypRedOdeSystemRhsAndJacobian<
    ind_var_type, reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType, HyperReducerType>;

  galerkin_system galSystem(trialSpace, fomSystem, hyperReducer);
  return ::pressio::ode::create_implicit_stepper(schemeName, std::move(galSystem));
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
&& RealValuedSemiDiscreteFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
&& std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
#endif
auto create_unsteady_implicit_problem(::pressio::ode::StepScheme schemeName,   /*(4)*/
				      const TrialSubspaceType & trialSpace,
				      const FomSystemType & fomSystem,
				      const MaskerType & masker,
				      const HyperReducerType & hyperReducer)
{

  impl::valid_scheme_for_implicit_galerkin_else_throw(schemeName, "galerkin_masked_implicit");

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using default_types      = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinMaskedOdeSystemRhsAndJacobian<
    ind_var_type, reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType,
    MaskerType, HyperReducerType>;

  galerkin_system galSystem(trialSpace, fomSystem, masker, hyperReducer);
  return ::pressio::ode::create_implicit_stepper(schemeName, std::move(galSystem));
}

// -------------------------------------------------------------
// fully-discrete
// -------------------------------------------------------------

template<
  std::size_t TotalNumberOfDesiredStates,
  class TrialSubspaceType,
  class FomSystemType>
#ifdef PRESSIO_ENABLE_CXX20
requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
&& RealValuedFullyDiscreteSystemWithJacobianAction<
    FomSystemType, TotalNumberOfDesiredStates, typename TrialSubspaceType::basis_matrix_type>
#endif
auto create_unsteady_implicit_problem(const TrialSubspaceType & trialSpace,    /*(5)*/
				      const FomSystemType & fomSystem)
{

  using ind_var_type = typename FomSystemType::time_type;
  using reduced_state_type    = typename TrialSubspaceType::reduced_state_type;
  using default_types         = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinDefaultFullyDiscreteSystem<
    TotalNumberOfDesiredStates, ind_var_type,
    reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType>;

  galerkin_system galSystem(trialSpace, fomSystem);
  return ::pressio::ode::create_implicit_stepper<
    TotalNumberOfDesiredStates>(std::move(galSystem));
}

template<
  std::size_t TotalNumberOfDesiredStates,
  class TrialSubspaceType,
  class FomSystemType,
  class HyperReducerType>
auto create_unsteady_implicit_problem(const TrialSubspaceType & trialSpace,
                                      const FomSystemType & fomSystem,
                                      const HyperReducerType & hyperReducer)
{

  using independent_variable_type = typename FomSystemType::time_type;
  using reduced_state_type        = typename TrialSubspaceType::reduced_state_type;
  using default_types             = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  // the "system" implements the math
  using galerkin_system = impl::GalerkinHypRedFullyDiscreteSystem<
    TotalNumberOfDesiredStates, independent_variable_type,
    reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType,
    HyperReducerType>;

  galerkin_system galSystem(trialSpace, fomSystem, hyperReducer);
  return ::pressio::ode::create_implicit_stepper<
    TotalNumberOfDesiredStates>(std::move(galSystem));
}

// -------------------------------------------------------------
// mixed fom/rom
// -------------------------------------------------------------
namespace experimental{

// because this is the mixed case FOM/ROM, we need to check that the FomSystem
// meets the API needed for ROM but also that needed for doing ode stepping directly on it
template<class TrialSubspaceType, class FomSystemType>
#ifdef PRESSIO_ENABLE_CXX20
requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
&& RealValuedSemiDiscreteFomWithJacobianAction<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
&& std::same_as<typename TrialSubspaceType::full_state_type, typename FomSystemType::state_type>
&& ::pressio::ode::RealValuedOdeSystemFusingRhsAndJacobian<FomSystemType>
#endif
auto create_unsteady_implicit_problem_mixed_fom(::pressio::ode::StepScheme schemeName,
						const TrialSubspaceType & trialSpace,
						const FomSystemType & fomSystem)
{

#if !defined PRESSIO_ENABLE_CXX20
  static_assert(PossiblyAffineTrialColumnSubspace<TrialSubspaceType>::value);
  static_assert(RealValuedSemiDiscreteFomWithJacobianAction<
		FomSystemType, typename TrialSubspaceType::basis_matrix_type>::value);
  static_assert(::pressio::ode::RealValuedOdeSystemFusingRhsAndJacobian<FomSystemType>::value);
#endif

  if (   schemeName != ::pressio::ode::StepScheme::BDF1
      && schemeName != ::pressio::ode::StepScheme::BDF2)
  {
    throw std::runtime_error("galerkin mixed fom/rom currently accepting BDF1 or BDF2");
  }

  using ind_var_type          = typename FomSystemType::time_type;
  using reduced_state_type    = typename TrialSubspaceType::reduced_state_type;
  using default_types         = ImplicitGalerkinDefaultReducedOperatorsTraits<reduced_state_type>;
  using reduced_residual_type = typename default_types::reduced_residual_type;
  using reduced_jacobian_type = typename default_types::reduced_jacobian_type;

  using galerkin_system = impl::GalerkinDefaultOdeSystemRhsAndJacobian<
    ind_var_type, reduced_state_type, reduced_residual_type,
    reduced_jacobian_type, TrialSubspaceType, FomSystemType>;

  using return_type = impl::GalerkinUnsteadyDefaultProblemRomFom<
    ind_var_type, FomSystemType, TrialSubspaceType, galerkin_system, impl::MixedFomRomStepper>;
  return return_type(schemeName, trialSpace, fomSystem);
}
}//end namespace experimental

}}} // end pressio::rom::galerkin
#endif  // ROM_GALERKIN_UNSTEADY_IMPLICIT_HPP_
