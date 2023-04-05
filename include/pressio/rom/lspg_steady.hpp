
#ifndef ROM_LSPG_STEADY_HPP_
#define ROM_LSPG_STEADY_HPP_

#include "./impl/lspg_steady_system_default.hpp"
#include "./impl/lspg_steady_system_masked.hpp"

namespace pressio{ namespace rom{ namespace lspg{

// -------------------------------------------------------------
// default or hyp-red
// impl-wise, default and hyp-red LSPG are the same
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType>
#ifdef PRESSIO_ENABLE_CXX20
  requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
  && RealValuedSteadyFomWithJacobianAction<
      FomSystemType, typename TrialSubspaceType::basis_matrix_type>
  && std::same_as<
      typename TrialSubspaceType::full_state_type,
      typename FomSystemType::state_type>
  //
  && (Traits<typename FomSystemType::state_type>::rank == 1)
  && (Traits<typename FomSystemType::residual_type>::rank == 1)
  && (Traits<typename TrialSubspaceType::basis_matrix_type>::rank == 2)
  && (Traits<impl::fom_jac_action_t<FomSystemType,
             typename TrialSubspaceType::basis_matrix_type> >::rank == 2)
#endif
auto create_steady_problem(const TrialSubspaceType & trialSpace,   /*(1)*/
			   const FomSystemType & fomSystem)

{

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = utils::NoOperation<void>;
  using system_type = impl::LspgSteadyDefaultSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, scaler_type>;
  return system_type(trialSpace, fomSystem, scaler_type{});
}

// -------------------------------------------------------------
// masked
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType>
#ifdef PRESSIO_ENABLE_CXX20
  requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
  && RealValuedSteadyFomWithJacobianAction<
      FomSystemType, typename TrialSubspaceType::basis_matrix_type>
  && std::same_as<
      typename TrialSubspaceType::full_state_type,
      typename FomSystemType::state_type>
  //
  && (Traits<typename FomSystemType::state_type>::rank == 1)
  && (Traits<typename FomSystemType::residual_type>::rank == 1)
  && (Traits<typename TrialSubspaceType::basis_matrix_type>::rank == 2)
  && (Traits<impl::fom_jac_action_t<FomSystemType,
             typename TrialSubspaceType::basis_matrix_type> >::rank == 2)
  //
  && MaskableWith<typename FomSystemType::residual_type, MaskerType>
  && MaskableWith<impl::fom_jac_action_on_trial_space_t<FomSystemType, TrialSubspaceType>, MaskerType>
#endif
auto create_steady_problem(const TrialSubspaceType & trialSpace,  /*(2)*/
			   const FomSystemType & fomSystem,
			   const MaskerType & masker)
{

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = utils::NoOperation<void>;
  using system_type = impl::LspgSteadyMaskedSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, MaskerType, scaler_type>;
  return system_type(trialSpace, fomSystem, masker, scaler_type{});
}


namespace experimental{
// -------------------------------------------------------------
// default or hyp-red with scaling
// impl-wise, default and hyp-red LSPG are the same
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class OperatorScalerType>
#ifdef PRESSIO_ENABLE_CXX20
  requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
  && RealValuedSteadyFomWithJacobianAction<
      FomSystemType, typename TrialSubspaceType::basis_matrix_type>
  && std::same_as<
      typename TrialSubspaceType::full_state_type,
      typename FomSystemType::state_type>
  //
  && (Traits<typename FomSystemType::state_type>::rank == 1)
  && (Traits<typename FomSystemType::residual_type>::rank == 1)
  && (Traits<typename TrialSubspaceType::basis_matrix_type>::rank == 2)
  && (Traits<impl::fom_jac_action_t<FomSystemType,
             typename TrialSubspaceType::basis_matrix_type> >::rank == 2)
  //
  && requires(const OperatorScalerType & scaler,
	      typename FomSystemType::state_type & s,
	      typename FomSystemType::residual_type & r,
	      std::optional<
	        impl::fom_jac_action_t<FomSystemType, typename TrialSubspaceType::basis_matrix_type> *
	      > jaO)
  {
    scaler(s, r, jaO);
  }
#endif
auto create_steady_problem(const TrialSubspaceType & trialSpace,
			   const FomSystemType & fomSystem,
			   const OperatorScalerType & scaler)
{

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = std::reference_wrapper<const OperatorScalerType>;
  using system_type = impl::LspgSteadyDefaultSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, scaler_type>;
  return system_type(trialSpace, fomSystem, scaler);
}


// -------------------------------------------------------------
// default or hyp-red with scaling
// impl-wise, default and hyp-red LSPG are the same
// -------------------------------------------------------------
template<
  class TrialSubspaceType,
  class FomSystemType,
  class MaskerType,
  class OperatorScalerType>
#ifdef PRESSIO_ENABLE_CXX20
  requires PossiblyAffineRealValuedTrialColumnSubspace<TrialSubspaceType>
  && RealValuedSteadyFomWithJacobianAction<
      FomSystemType, typename TrialSubspaceType::basis_matrix_type>
  && std::same_as<
      typename TrialSubspaceType::full_state_type,
      typename FomSystemType::state_type>
  //
  && (Traits<typename FomSystemType::state_type>::rank == 1)
  && (Traits<typename FomSystemType::residual_type>::rank == 1)
  && (Traits<typename TrialSubspaceType::basis_matrix_type>::rank == 2)
  && (Traits<impl::fom_jac_action_t<FomSystemType,
             typename TrialSubspaceType::basis_matrix_type> >::rank == 2)
  //
  && MaskableWith<typename FomSystemType::residual_type, MaskerType>
  && MaskableWith<
      impl::fom_jac_action_on_trial_space_t<FomSystemType, TrialSubspaceType>, MaskerType>
  //
  && requires(const OperatorScalerType & scaler,
	      const typename FomSystemType::state_type & s,
	      impl::mask_action_t<MaskerType, typename FomSystemType::residual_type> & maskedR,
	      std::optional<
	       impl::mask_action_t<
	        MaskerType,
	        impl::fom_jac_action_t<FomSystemType, typename TrialSubspaceType::basis_matrix_type>
	       > *> maskedJaO)
  {
    scaler(s, maskedR, maskedJaO);
  }
#endif
auto create_steady_problem(const TrialSubspaceType & trialSpace,
			   const FomSystemType & fomSystem,
			   const MaskerType & masker,
			   const OperatorScalerType & scaler)
{

  using reduced_state_type = typename TrialSubspaceType::reduced_state_type;
  using scaler_type = std::reference_wrapper<const OperatorScalerType>;
  using system_type = impl::LspgSteadyMaskedSystem<
    reduced_state_type, TrialSubspaceType, FomSystemType, MaskerType, scaler_type>;
  return system_type(trialSpace, fomSystem, masker, scaler);
}

} // end experimental

}}} // end pressio::rom::lspg
#endif  // ROM_LSPG_STEADY_HPP_
