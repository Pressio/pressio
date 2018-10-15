
#ifndef ROM_LSPG_RESIDUAL_POLICY_HPP_
#define ROM_LSPG_RESIDUAL_POLICY_HPP_

#include "rom_forward_declarations.hpp"
#include "../../../CORE_ALL"
#include "../../../ode/src/implicit/ode_residual_impl.hpp"
#include "../../../ode/src/implicit/policies/base/ode_implicit_residual_policy_base.hpp"
#include "rom_incremental_solution_base.hpp"

namespace rompp{ namespace rom{ namespace exp{

template<typename app_state_w_type,
	 typename app_res_w_type,
	 // basis operator type
	 typename phi_op_type,  
	 // weighting matrix by default is = void
	 typename A_type>
class RomLSPGResidualPolicy<::rompp::ode::ImplicitSteppersEnum::Euler,
			     app_state_w_type, app_res_w_type,
			     phi_op_type, A_type>
  : public ode::policy::ImplicitResidualPolicyBase<
                RomLSPGResidualPolicy<::rompp::ode::ImplicitSteppersEnum::Euler,
				       app_state_w_type, app_res_w_type,
				       phi_op_type, A_type>, 1, 0>,
    private IncrementalSolutionBase<
                RomLSPGResidualPolicy<::rompp::ode::ImplicitSteppersEnum::Euler,
				       app_state_w_type, app_res_w_type,
				       phi_op_type, A_type>, app_state_w_type>{
  
  using this_t = RomLSPGResidualPolicy<
    ::rompp::ode::ImplicitSteppersEnum::Euler,
    app_state_w_type, app_res_w_type, phi_op_type, A_type>;

  using base_pol_t = ::rompp::ode::policy::ImplicitResidualPolicyBase<this_t, 1, 0>;
  using base_incr_sol_t = rompp::rom::exp::IncrementalSolutionBase<this_t, app_state_w_type>;
  using scalar_type = typename core::details::traits<app_state_w_type>::scalar_t;
  
 private:
  mutable app_res_w_type appRHS_;
  phi_op_type * phi_ = nullptr;
  A_type * A_ = nullptr;
  
  using base_incr_sol_t::y0FOM_;
  using base_incr_sol_t::yFOM_;
  
public:
  template <typename T=A_type,
   core::meta::enable_if_t<std::is_void<T>::value> * = nullptr>
  RomLSPGResidualPolicy(const app_state_w_type & y0fom,
  			const app_res_w_type & r0fom,
  			phi_op_type & phiOp)
    : base_incr_sol_t(y0fom), appRHS_(r0fom),
    phi_(&phiOp), yFOMnm1_(y0fom){}

  template <typename T=A_type,
   core::meta::enable_if_t<std::is_void<T>::value> * = nullptr>
  RomLSPGResidualPolicy(const app_state_w_type & y0fom,
  			phi_op_type & phiOp)
    : base_incr_sol_t(y0fom), phi_(&phiOp), yFOMnm1_(y0fom){}
  
  RomLSPGResidualPolicy() = delete;
  ~RomLSPGResidualPolicy() = default;
  //----------------------------------------------------------------
  
  template <typename app_t>
  void operator()(const app_state_w_type & odeY,
  		  app_res_w_type & odeR,
  		  const std::array<app_state_w_type, 1> & oldYs,
  		  const app_t & app,
  		  scalar_type t,
  		  scalar_type dt) const{
    // here y_n = odeY is reduced state
    // need to compute: R( phi y_n) = phi y_n - phi y_n-1 - dt * f(phi y)
    
    // odeY is the REDUCED state, we need to reconstruct FOM state
    phi_->apply(odeY, yFOM_);

    // reconstruct FOM state at previous step n-1
    // the previous state is stored inside oldYs
    phi_->apply(oldYs[0], yFOMnm1_);
    
    // since we are advancing the Incremental Solution,
    // to compute the app residual we need to add the
    // FOM initial condition to get full state 
    yFOM_ += (*y0FOM_);
    yFOMnm1_ += (*y0FOM_);
    
    /// query the application for the SPACE residual 
    app.residual(*yFOM_.data(), *odeR.data(), t);

    // do time discrete residual
    ode::impl::implicit_euler_time_discrete_residual(yFOM_, yFOMnm1_, odeR, dt);
    
    // /// apply weighting
    // if (A_)
    //   A_->applyTranspose(appRHS_, odeR);
  }
  //----------------------------------------------------------------

  template <typename app_t>
  app_res_w_type operator()(const app_state_w_type & odeY,
			    const std::array<app_state_w_type, 1> & oldYs,
			    const app_t & app,
			    scalar_type t,
			    scalar_type dt) const{
    // here y_n = odeY is reduced state
    // need to compute: R( phi y_n) = phi y_n - phi y_n-1 - dt * f(phi y)
    
    // odeY is the REDUCED state, we need to reconstruct FOM state
    phi_->apply(odeY, yFOM_);

    // reconstruct FOM state at previous step n-1
    // the previous state is stored inside oldYs
    phi_->apply(oldYs[0], yFOMnm1_);
    
    // since we are advancing the Incremental Solution,
    // to compute the app residual we need to add the
    // FOM initial condition to get full state 
    yFOM_ += (*y0FOM_);
    yFOMnm1_ += (*y0FOM_);
    
    /// query the application for the SPACE residual 
    auto RR = app.residual(*yFOM_.data(), t);
    app_res_w_type RRw(RR);

    // do time discrete residual
    ode::impl::implicit_euler_time_discrete_residual(yFOM_, yFOMnm1_,
						     RRw, dt);
    /// apply weighting
    // if (A_)
    //   A_->applyTranspose(appRHS_, odeR);
    
    return RRw;
  }
  //----------------------------------------------------------------

private:
  friend base_pol_t;
  friend base_incr_sol_t;
  mutable app_state_w_type yFOMnm1_;

};//end class
  
}}}//end namespace rompp::rom::exp
#endif 
