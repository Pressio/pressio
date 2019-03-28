
#ifndef ROM_LSPG_STEADY_JACOBIAN_POLICY_HPP_
#define ROM_LSPG_STEADY_JACOBIAN_POLICY_HPP_

#include "../../rom_forward_declarations.hpp"
#include "../../rom_data_fom_states.hpp"

namespace rompp{ namespace rom{

template<
  typename fom_states_data,
  typename apply_jac_return_type,
  typename fom_apply_jac_policy,
  typename decoder_type
  >
class LSPGSteadyJacobianPolicy
  : protected fom_states_data,
    protected fom_apply_jac_policy{

protected:
  using this_t = LSPGSteadyJacobianPolicy<
  fom_states_data, apply_jac_return_type,
  fom_apply_jac_policy, decoder_type>;

  using fom_states_data::yFom_;

public:
  static constexpr bool isResidualPolicy_ = false;
  using apply_jac_return_t = apply_jac_return_type;

public:
  LSPGSteadyJacobianPolicy() = delete;

  ~LSPGSteadyJacobianPolicy() = default;

  LSPGSteadyJacobianPolicy(const fom_states_data	& fomStates,
			   const fom_apply_jac_policy	& applyJacFunctor,
			   const apply_jac_return_type	& applyJacObj,
			   const decoder_type		& decoder)
    : fom_states_data(fomStates),
      fom_apply_jac_policy(applyJacFunctor),
      JJ_(applyJacObj),
      decoderObj_(decoder){}

public:

  template <typename lspg_state_t,
	    typename lspg_jac_t,
	    typename app_t>
  void operator()(const lspg_state_t & romY,
		  lspg_jac_t	     & romJJ,
  		  const app_t	     & app) const
  {
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("lspg apply jac");
#endif

    // todo: this is not needed if jacobian is called after resiudal
    // because residual takes care of reconstructing the fom state
    //    timer->start("reconstruct fom state");
    fom_states_data::template reconstructCurrentFomState(romY);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("fom apply jac");
    const auto & basis = decoderObj_.getReferenceToJacobian();
    fom_apply_jac_policy::evaluate(app, yFom_, basis, romJJ);
    timer->stop("fom apply jac");
#else
    const auto & basis = decoderObj_.getReferenceToJacobian();
    fom_apply_jac_policy::evaluate(app, yFom_, basis, romJJ);
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("lspg apply jac");
#endif
  }

  template <typename lspg_state_t, typename app_t>
  apply_jac_return_t operator()(const lspg_state_t & romY,
				const app_t	   & app) const
  {
    (*this).template operator()(romY, JJ_, app);
    return JJ_;
  }

protected:
  mutable apply_jac_return_t JJ_   = {};
  const decoder_type & decoderObj_ = {};

};

}}//end namespace rompp::rom
#endif
