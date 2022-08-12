
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
  class TrialSpaceType,
  class FomSystemType
  >
class LspgSteadyDefaultSystem
{

  // need to deduce the type of the action of the fom jacobian
  // which becomes the jacobian_type of the problem
  using fom_jac_action_result_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
	     (std::declval<typename TrialSpaceType::basis_type const &>())
	     );

public:
  // aliases required by the pressio solvers
  using state_type    = ReducedStateType;
  using residual_type = typename FomSystemType::residual_type;
  using jacobian_type = fom_jac_action_result_type;

  LspgSteadyDefaultSystem() = delete;

  LspgSteadyDefaultSystem(const TrialSpaceType & trialSpace,
			  const FomSystemType & fomSystem)
    : trialSpace_(trialSpace), fomSystem_(fomSystem),
      fomState_(fomSystem.createState())
  {}

public:
  jacobian_type createJacobian() const{
    return fomSystem_.get().createApplyJacobianResult(trialSpace_.get().viewBasis());
  }

  residual_type createResidual() const{
    return fomSystem_.get().createResidual();
  }

  void residualAndJacobian(const state_type & lspgState,
			   residual_type & R,
			   jacobian_type & J,
			   bool recomputeJacobian) const
  {
    trialSpace_.get().mapFromReducedState(lspgState, fomState_);
    fomSystem_.get().residual(fomState_, R);

    if (recomputeJacobian){
      const auto & phi = trialSpace_.get().viewBasis();
      fomSystem_.get().applyJacobian(fomState_, phi, J);
    }
  }

private:
  std::reference_wrapper<const TrialSpaceType> trialSpace_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
};


}}} // end pressio::rom::impl
#endif
