
#ifndef ROM_LSPG_RESIDUAL_POLICY_HPP_
#define ROM_LSPG_RESIDUAL_POLICY_HPP_

#include "../rom_forward_declarations.hpp"
#include "../../../CORE_ALL"
#include "../../../ode/src/implicit/ode_residual_impl.hpp"
#include "../../../ode/src/implicit/policies/base/ode_implicit_residual_policy_base.hpp"
#include "../rom_data_base.hpp"

namespace rompp{ namespace rom{

template< typename app_state_w_type,
	  typename app_res_w_type,
	  typename phi_op_type,
	  int maxNstates,
	  int maxNrhs,
	  typename A_type>
class RomLSPGResidualPolicy
  : public ode::policy::ImplicitResidualPolicyBase<
                RomLSPGResidualPolicy<app_state_w_type,
				      app_res_w_type,
				      phi_op_type,
				      maxNstates, maxNrhs,
				      A_type>>,
    protected RomStateData<app_state_w_type, phi_op_type,  maxNstates>,
    protected RomRHSData<app_res_w_type, maxNrhs>
{

  using this_t 		= RomLSPGResidualPolicy<app_state_w_type,
						app_res_w_type,
			                        phi_op_type,
						maxNstates, maxNrhs,
						A_type>;

  using base_pol_t	  = ode::policy::ImplicitResidualPolicyBase<this_t>;
  using base_state_data_t = rom::RomStateData<app_state_w_type, phi_op_type, maxNstates>;
  using base_rhs_data_t   = rom::RomRHSData<app_res_w_type, maxNrhs>;
  friend base_pol_t;

public:
  using app_res_w_t = app_res_w_type;
    using scalar_type	  = typename core::details::traits<app_state_w_type>::scalar_t;

 private:
  using base_state_data_t::y0FOM_;
  using base_state_data_t::yFOM_;
  using base_state_data_t::yFOMold_;
  using base_state_data_t::phi_;
  using base_rhs_data_t::appRHS_;
  //  using base_data_t::A_;

public:

  RomLSPGResidualPolicy() = delete;

  ~RomLSPGResidualPolicy() = default;

  template <typename T=A_type,
	    core::meta::enable_if_t<std::is_void<T>::value> * = nullptr>
  RomLSPGResidualPolicy(const app_state_w_type & y0fom,
			const app_res_w_type & r0fom,
			phi_op_type & phiOp)
    : base_state_data_t(y0fom, phiOp), base_rhs_data_t(r0fom){}

 public:
  template <::rompp::ode::ImplicitEnum odeMethod,
	     int numAuxStates,
	     typename ode_state_t,
	     typename app_t>
  app_res_w_type operator()(const ode_state_t & odeY,
			    const std::array<ode_state_t,numAuxStates> & oldYs,
			    const app_t & app,
			    scalar_type t,
			    scalar_type dt) const
  {

#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("lspg reconstruct FOM state");
#endif
    // odeY, oldYsa are REDUCED states, reconstruct FOM states
    base_state_data_t::template reconstructFOMStates<numAuxStates>(odeY, oldYs);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("lspg reconstruct FOM state");
#endif


#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("lspg space residual");
#endif
    /// query the application for the SPACE residual
    app.residual(*yFOM_.data(), *appRHS_[0].data(), t);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("lspg space residual");
#endif


#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("lspg time discrete residual");
#endif
    /// do time discrete residual
    ode::impl::implicit_time_discrete_residual<odeMethod, maxNstates>(yFOM_,
							     yFOMold_,
							     appRHS_[0],
							     dt);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("lspg time discrete residual");
#endif

    // /// apply weighting:
    // if (A_) A_->applyTranspose(appRHS_, odeR);
    return appRHS_[0];
  }
  //----------------------------------------------------------------


  template <::rompp::ode::ImplicitEnum odeMethod,
	     int numAuxStates,
	     typename ode_state_t,
	     typename ode_res_t,
	     typename app_t>
  void operator()(const ode_state_t & odeY,
  		  ode_res_t & odeR,
  		  const std::array<ode_state_t,numAuxStates> & oldYs,
  		  const app_t & app,
  		  scalar_type t,
  		  scalar_type dt) const
  {
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("lspg reconstruct FOM state");
#endif
    // odeY, oldYsa are REDUCED states, reconstruct FOM states
    base_state_data_t::template reconstructFOMStates<numAuxStates>(odeY, oldYs);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("lspg reconstruct FOM state");
#endif


#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("lspg space residual");
#endif
    app.residual(*yFOM_.data(), *odeR.data(), t);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("lspg space residual");
#endif


#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("lspg time discrete residual");
#endif
    ode::impl::implicit_time_discrete_residual<odeMethod, maxNstates>(yFOM_,
								      yFOMold_,
								      odeR,
								      dt);
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("lspg time discrete residual");
#endif

    // /// apply weighting
    // if (A_) A_->applyTranspose(appRHS_, odeR);
  }


};//end class

}}//end namespace rompp::rom
#endif
