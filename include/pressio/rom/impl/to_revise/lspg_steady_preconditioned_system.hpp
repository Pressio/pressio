
#ifndef ROM_IMPL_TO_REVISE_LSPG_STEADY_PRECONDITIONED_SYSTEM_HPP_
#define ROM_IMPL_TO_REVISE_LSPG_STEADY_PRECONDITIONED_SYSTEM_HPP_

namespace pressio{ namespace rom{ namespace impl{

template <
  class LspgStateType, class TrialSpaceType,
  class FomSystemType, class PreconditionerType
  >
class LspgSteadyPreconditionedSystem
  : public LspgSteadyDefaultSystem<
  LspgStateType, TrialSpaceType, FomSystemType>

{
  using base_t = LspgSteadyDefaultSystem<
    LspgStateType, TrialSpaceType, FomSystemType>;

public:
  // required aliases
  using typename base_t::state_type;
  using typename base_t::residual_type;
  using typename base_t::jacobian_type;

  LspgSteadyPreconditionedSystem() = delete;

  LspgSteadyPreconditionedSystem(TrialSpaceType & space,
         const FomSystemType & fomSystem,
         const PreconditionerType & preconditioner)
    : base_t(space, fomSystem), prec_(preconditioner){}

public:
  jacobian_type createJacobian() const{
    return base_t::createJacobian();
  }

  residual_type createResidual() const{
    return base_t::createResidual();
  }

  void residualAndJacobian(const state_type & lspgState,
        residual_type & R,
        jacobian_type & J,
        bool recomputeJacobian = true) const
  {
    base_t::residualAndJacobian(lspgState, R, J, recomputeJacobian);
    prec_(base_t::fomState_, R);
    if (recomputeJacobian){
      prec_(base_t::fomState_, J);
    }
  }

private:
  std::reference_wrapper<const PreconditionerType> prec_;
};

}}} // end pressio::rom::impl
#endif  // ROM_IMPL_TO_REVISE_LSPG_STEADY_PRECONDITIONED_SYSTEM_HPP_
