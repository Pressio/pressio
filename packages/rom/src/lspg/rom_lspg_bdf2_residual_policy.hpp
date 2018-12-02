
#ifndef ROM_LSPG_BDF2_RESIDUAL_POLICY_HPP_
#define ROM_LSPG_BDF2_RESIDUAL_POLICY_HPP_

#include "../rom_forward_declarations.hpp"
#include "../../../CORE_ALL"
#include "../../../ode/src/implicit/ode_residual_impl.hpp"
#include "../../../ode/src/implicit/policies/base/ode_implicit_residual_policy_base.hpp"
#include "../rom_incremental_solution_base.hpp"

namespace rompp{ namespace rom{

template<typename app_state_w_type,
	 typename app_res_w_type,
	 typename phi_op_type,
	 typename A_type>
class RomLSPGResidualPolicy<::rompp::ode::ImplicitSteppersEnum::BDF2,
			     app_state_w_type, app_res_w_type,
			     phi_op_type, A_type>
  : public ode::policy::ImplicitResidualPolicyBase<
                RomLSPGResidualPolicy<::rompp::ode::ImplicitSteppersEnum::BDF2,
				       app_state_w_type, app_res_w_type,
				       phi_op_type, A_type>, 2, 0>,
    private IncrementalSolutionBase<
		RomLSPGResidualPolicy<::rompp::ode::ImplicitSteppersEnum::BDF2,
				       app_state_w_type, app_res_w_type,
				       phi_op_type, A_type>, app_state_w_type>{

  using this_t 		= RomLSPGResidualPolicy<::rompp::ode::ImplicitSteppersEnum::BDF2,
    			  app_state_w_type, app_res_w_type, phi_op_type, A_type>;
  using base_pol_t 	= ::rompp::ode::policy::ImplicitResidualPolicyBase<this_t, 2, 0>;
  using base_incr_sol_t = rompp::rom::IncrementalSolutionBase<this_t, app_state_w_type>;
  using scalar_type 	= typename core::details::traits<app_state_w_type>::scalar_t;

 private:
  mutable app_res_w_type appRHS_ = {};
  phi_op_type * phi_ 		 = nullptr;
  A_type * A_ 			 = nullptr;

  using base_incr_sol_t::y0FOM_;
  using base_incr_sol_t::yFOM_;

 public:
  template <typename T=A_type,
	    core::meta::enable_if_t<std::is_void<T>::value> * = nullptr>
    RomLSPGResidualPolicy(const app_state_w_type & y0fom,
			  const app_res_w_type & r0fom,
			  phi_op_type & phiOp)
    : base_incr_sol_t(y0fom), appRHS_(r0fom),
    phi_(&phiOp), yFOMnm1_(y0fom), yFOMnm2_(y0fom){}

  template <typename T=A_type,
   core::meta::enable_if_t<std::is_void<T>::value> * = nullptr>
  RomLSPGResidualPolicy(const app_state_w_type & y0fom,
  			phi_op_type & phiOp)
    : base_incr_sol_t(y0fom), phi_(&phiOp),
    yFOMnm1_(y0fom), yFOMnm2_(y0fom){}

  RomLSPGResidualPolicy() = delete;

  ~RomLSPGResidualPolicy() = default;

  //----------------------------------------------------------------

  template <typename ode_state_t,
	    typename app_t>
  app_res_w_type operator()(const ode_state_t & odeY,
			    const std::array<ode_state_t, 2> & oldYs,
			    const app_t & app,
			    scalar_type t,
			    scalar_type dt) const
  {
    // odeY, oldYsa are REDUCED states, reconstruct FOM states
    reconstructFOMStates(odeY, oldYs[1], oldYs[0]);

    /// query the application for the SPACE residual
    app.residual(*yFOM_.data(), *appRHS_.data(), t);

    /// do time discrete residual
    ode::impl::implicit_bdf2_time_discrete_residual(yFOM_, yFOMnm1_,
						    yFOMnm2_, appRHS_, dt);
    // /// apply weighting:
    // if (A_) A_->applyTranspose(appRHS_, odeR);
    return appRHS_;
  }
  //----------------------------------------------------------------


  template <typename ode_state_t,
	    typename ode_res_t,
	    typename app_t>
  void operator()(const ode_state_t & odeY,
  		  ode_res_t & odeR,
  		  const std::array<ode_state_t, 2> & oldYs,
  		  const app_t & app,
  		  scalar_type t,
  		  scalar_type dt) const
  {
    reconstructFOMStates(odeY, oldYs[1], oldYs[0]);
    app.residual(*yFOM_.data(), *odeR.data(), t);
    ode::impl::implicit_bdf2_time_discrete_residual(yFOM_, yFOMnm1_,
						    yFOMnm2_, odeR, dt);
    // /// apply weighting
    // if (A_) A_->applyTranspose(appRHS_, odeR);
  }
  //----------------------------------------------------------------

 private:
  template <typename ode_state_t>
  void reconstructFOMStates(const ode_state_t & odeY,
 			    const ode_state_t & odeYm1,
			    const ode_state_t & odeYm2) const {
    phi_->apply(odeY, yFOM_);
    phi_->apply(odeYm1, yFOMnm1_);
    phi_->apply(odeYm2, yFOMnm2_);
    yFOM_ += (*y0FOM_);
    yFOMnm1_ += (*y0FOM_);
    yFOMnm2_ += (*y0FOM_);
  }
  //----------------------------------------------------------------

 private:
  friend base_pol_t;
  friend base_incr_sol_t;
  mutable app_state_w_type yFOMnm1_ = {};
  mutable app_state_w_type yFOMnm2_ = {};

};//end class

}}//end namespace rompp::rom
#endif
