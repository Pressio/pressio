
#ifndef ROM_LSPG_JACOBIAN_POLICY_HPP_
#define ROM_LSPG_JACOBIAN_POLICY_HPP_

#include "../../rom_forward_declarations.hpp"
#include "rom_lspg_time_discrete_jacobian.hpp"
#include "../../../../ode/src/implicit/policies/base/ode_jacobian_policy_base.hpp"
#include "../../rom_data_fom_states.hpp"

namespace rompp{ namespace rom{

template<
  typename fom_states_data,
  typename apply_jac_return_type,
  typename fom_apply_jac_policy,
  typename decoder_type,
  typename td_ud_ops
  >
class LSPGJacobianPolicy
  : public ode::policy::JacobianPolicyBase<
	LSPGJacobianPolicy<fom_states_data,
			   apply_jac_return_type,
			   fom_apply_jac_policy,
			   decoder_type,
			   td_ud_ops>>,
    protected fom_states_data,
    protected fom_apply_jac_policy{

protected:
  using this_t = LSPGJacobianPolicy<fom_states_data,
				    apply_jac_return_type,
				    fom_apply_jac_policy,
				    decoder_type,
				    td_ud_ops>;

  friend ode::policy::JacobianPolicyBase<this_t>;
  using fom_states_data::yFom_;
  using fom_states_data::yFomOld_;

  const td_ud_ops & tdOps_;

public:
  static constexpr bool isResidualPolicy_ = false;
  using apply_jac_return_t = apply_jac_return_type;

public:
  LSPGJacobianPolicy() = delete;
  ~LSPGJacobianPolicy() = default;

  LSPGJacobianPolicy(const fom_states_data	 & fomStates,
		     const fom_apply_jac_policy  & applyJacFunctor,
		     const apply_jac_return_type & applyJacObj,
		     const decoder_type		 & decoder)
    : fom_states_data(fomStates),
      fom_apply_jac_policy(applyJacFunctor),
      JJ_(applyJacObj),
      decoderObj_(decoder){}


  LSPGJacobianPolicy(const fom_states_data	 & fomStates,
		     const fom_apply_jac_policy  & applyJacFunctor,
		     const apply_jac_return_type & applyJacObj,
		     const decoder_type		 & decoder,
		     const td_ud_ops & tdOps)
    : fom_states_data(fomStates),
      fom_apply_jac_policy(applyJacFunctor),
      JJ_(applyJacObj),
      decoderObj_(decoder),
      tdOps_(tdOps){}

  // user-defined ops is void
  template <
    ode::ImplicitEnum odeMethod,
    typename T = td_ud_ops,
    rompp::mpl::enable_if_t<
      std::is_void<T>::value == true
      > * = nullptr
    >
  struct td_dispatcher{
    template <typename ... Args>
    static void evaluate(Args && ... args){
      rom::impl::time_discrete_jacobian
	<odeMethod>(std::forward<Args>(args)...);
    }
  };

  // user defined ops is non-trivial
  template <
    ode::ImplicitEnum odeMethod,
    typename T = td_ud_ops,
    rompp::mpl::enable_if_t<
      std::is_void<T>::value == false
      > * = nullptr
    >
  struct td_dispatcher{
    template <typename ... Args>
    static void evaluate(Args && ... args){
      rom::impl::time_discrete_jacobian
	<odeMethod>(tdOps_, std::forward<Args>(args)...);
    }
  };

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
#endif
    const auto & basis = decoderObj_.getReferenceToJacobian();
    fom_apply_jac_policy::evaluate(app, yFom_, basis, romJJ, t);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("fom apply jac");
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("time discrete jacob");
#endif
    this->td_dispatcher<odeMethod>(romJJ, dt, basis);
    //rom::impl::time_discrete_jacobian<odeMethod>(romJJ, dt, basis);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("time discrete jacob");
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
  mutable apply_jac_return_t JJ_	= {};
  const decoder_type & decoderObj_	= {};

};

}}//end namespace rompp::rom
#endif
