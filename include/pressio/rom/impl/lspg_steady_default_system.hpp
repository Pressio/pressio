
#ifndef PRESSIO_ROM_IMPL_LSPG_STEADY_DEFAULT_SYSTEM_HPP_
#define PRESSIO_ROM_IMPL_LSPG_STEADY_DEFAULT_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <class LspgStateType, class TrialSpaceType, class FomSystemType>
class LspgSteadyDefaultSystem
{

public:
  // required aliases
  using state_type    = LspgStateType;
  using residual_type = typename FomSystemType::residual_type;
  // need to deduce the type of the action of the fom jacobian
  using jacobian_type =
    decltype(std::declval<FomSystemType const>().createApplyJacobianResult
      (std::declval<typename TrialSpaceType::basis_type const &>())
      );

  LspgSteadyDefaultSystem() = delete;

  LspgSteadyDefaultSystem(TrialSpaceType & space,
       const FomSystemType & fomSystem)
    : space_(space), fomSystem_(fomSystem),
      fomState_(fomSystem.createState())
  {}

public:
  jacobian_type createJacobian() const{
    return fomSystem_.get().createApplyJacobianResult(space_.get().viewBasis());
  }

  residual_type createResidual() const{
    return fomSystem_.get().createResidual();
  }

  void residualAndJacobian(const state_type & lspgState,
        residual_type & R,
        jacobian_type & J,
        bool recomputeJacobian = true) const
  {
    space_.get().mapFromReducedState(lspgState, fomState_);
    fomSystem_.get().residual(fomState_, R);

    if (recomputeJacobian){
      const auto & phi = space_.get().viewBasis();
      fomSystem_.get().applyJacobian(fomState_, phi, J);
    }
  }

protected:
  std::reference_wrapper<TrialSpaceType> space_;
  std::reference_wrapper<const FomSystemType> fomSystem_;
  mutable typename FomSystemType::state_type fomState_;
};


}}} // end pressio::rom::impl
#endif
