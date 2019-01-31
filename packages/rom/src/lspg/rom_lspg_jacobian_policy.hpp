
#ifndef ROM_LSPG_JACOBIAN_POLICY_HPP_
#define ROM_LSPG_JACOBIAN_POLICY_HPP_

#include "../rom_forward_declarations.hpp"
#include "rom_lspg_time_discrete_jacobian.hpp"
#include "../../../ode/src/implicit/policies/base/ode_jacobian_policy_base.hpp"
#include "../rom_data_fom_states.hpp"

namespace rompp{ namespace rom{

template< typename fom_states_data,
	  typename apply_jac_return_type,
	  typename fom_apply_jac_policy>
class LSPGJacobianPolicy
  : public ode::policy::JacobianPolicyBase<
	LSPGJacobianPolicy<fom_states_data,
			   apply_jac_return_type,
			   fom_apply_jac_policy>>,
    protected fom_states_data,
    protected fom_apply_jac_policy{

protected:
  using this_t = LSPGJacobianPolicy<fom_states_data,
				    apply_jac_return_type,
				    fom_apply_jac_policy>;
  friend ode::policy::JacobianPolicyBase<this_t>;

  using fom_states_data::yFom_;
  using fom_states_data::yFomOld_;
  using fom_states_data::maxNstates_;
  using fom_states_data::decoderObj_;

public:
  static constexpr bool isResidualPolicy_ = false;
  using apply_jac_return_t = apply_jac_return_type;

public:
  LSPGJacobianPolicy() = delete;
  ~LSPGJacobianPolicy() = default;
  LSPGJacobianPolicy(const fom_states_data & fomStates,
		     const fom_apply_jac_policy & applyJacFunctor,
		     const apply_jac_return_type & applyJacObj)
    : fom_states_data(fomStates),
      fom_apply_jac_policy(applyJacFunctor),
      JJ_(applyJacObj){}

public:
  template <ode::ImplicitEnum odeMethod,
	    typename lspg_state_t,
	    typename lspg_jac_t,
	    typename app_t,
	    typename scalar_t>
  void operator()(const lspg_state_t & romY,
		  lspg_jac_t	     & romJJ,
  		  const app_t	     & app,
		  scalar_t	     t,
		  scalar_t	     dt) const
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
    const auto & basis = decoderObj_.getJacobianRef();
    fom_apply_jac_policy::evaluate(app, yFom_, basis, romJJ, t);
    timer->stop("fom apply jac");
#else
    const auto & basis = decoderObj_.getJacobianRef();
    fom_apply_jac_policy::evaluate(app, yFom_, basis, romJJ, t);
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("time discrete jacob");
    rom::impl::time_discrete_jacobian<odeMethod>(romJJ, dt, basis);
    timer->stop("time discrete jacob");
#else
    rom::impl::time_discrete_jacobian<odeMethod>(romJJ, dt, basis);
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("lspg apply jac");
#endif
  }

  template <ode::ImplicitEnum odeMethod,
	    typename lspg_state_t,
	    typename app_t,
	    typename scalar_t>
  apply_jac_return_t operator()(const lspg_state_t & romY,
				const app_t	   & app,
				scalar_t	   t,
				scalar_t	   dt) const
  {
    (*this).template operator()<odeMethod>(romY, JJ_, app, t, dt);
    return JJ_;
  }

protected:
  mutable apply_jac_return_t JJ_ = {};
};

}}//end namespace rompp::rom
#endif
