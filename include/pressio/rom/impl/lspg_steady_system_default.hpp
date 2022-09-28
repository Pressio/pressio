
#ifndef PRESSIO_ROM_IMPL_LSPG_STEADY_DEFAULT_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_LSPG_STEADY_DEFAULT_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

/*
LSPG steady default represents:

  min_x ||fom_r(phi x)||

- fom_r is the fom "residual"
- phi is the basis
*/
template <
  class ReducedStateType,
  class TrialSubspaceType,
  class FomSystemType
  >
class LspgSteadyDefaultSystem
{

  // need to deduce the type of the action of the fom jacobian
  // which becomes the jacobian_type of the problem
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<typename TrialSubspaceType::basis_matrix_type const &>())
	     );

public:
  // aliases required by the pressio solvers
  using state_type    = ReducedStateType;
  using residual_type = typename FomSystemType::residual_type;
  using jacobian_type = fom_jac_action_result_type;

  LspgSteadyDefaultSystem() = delete;

  LspgSteadyDefaultSystem(const TrialSubspaceType & trialSubspace,
			  const FomSystemType & fomSystem)
    : trialSubspace_(trialSubspace),
      fomSystem_(fomSystem),
      fomState_(trialSubspace.createFullState())
  {}

public:
  state_type createState() const{
    return trialSubspace_.get().createReducedState();
  }

  residual_type createResidual() const{
    return fomSystem_.get().createResidual();
  }

  jacobian_type createJacobian() const{
    return fomSystem_.get().createApplyJacobianResult(trialSubspace_.get().basisOfTranslatedSpace());
  }

  void residualAndJacobian(const state_type & lspgState,
			   residual_type & lsgpResidual,
			   jacobian_type & lspgJacobian,
			   bool computeJacobian) const
  {
    trialSubspace_.get().mapFromReducedState(lspgState, fomState_);
    fomSystem_.get().residual(fomState_, lsgpResidual);

    if (computeJacobian){
      const auto & phi = trialSubspace_.get().basisOfTranslatedSpace();
      fomSystem_.get().applyJacobian(fomState_, phi, lspgJacobian);
    }
  }

private:
  std::reference_wrapper<const TrialSubspaceType> trialSubspace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
};


}}} // end pressio::rom::impl
#endif
