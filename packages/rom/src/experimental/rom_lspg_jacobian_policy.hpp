
#ifndef ROM_LSPG_JACOBIAN_POLICY_HPP_
#define ROM_LSPG_JACOBIAN_POLICY_HPP_

#include "rom_forward_declarations.hpp"
#include "../../../CORE_ALL"
#include "../../../ode/src/implicit/ode_jacobian_impl.hpp"
#include "../../../ode/src/implicit/policies/base/ode_jacobian_policy_base.hpp"
#include "rom_incremental_solution_base.hpp"

namespace rompp{ namespace rom{ namespace exp{

template<typename app_state_w_type,
	 typename app_jac_w_type,
	 // basis operator type
	 typename phi_op_type,  
	 // weighting matrix by default is = void
	 typename A_type>
class RomLSPGJacobianPolicy<::rompp::ode::ImplicitSteppersEnum::Euler,
			     app_state_w_type, app_jac_w_type,
			     phi_op_type, A_type>
  : public ode::policy::JacobianPolicyBase<
                RomLSPGJacobianPolicy<::rompp::ode::ImplicitSteppersEnum::Euler,
				       app_state_w_type, app_jac_w_type,
				       phi_op_type, A_type>>,
    private IncrementalSolutionBase<
                RomLSPGJacobianPolicy<::rompp::ode::ImplicitSteppersEnum::Euler,
				       app_state_w_type, app_jac_w_type,
				       phi_op_type, A_type>, app_state_w_type>{
  
  using this_t = RomLSPGJacobianPolicy<
    ::rompp::ode::ImplicitSteppersEnum::Euler,
     app_state_w_type, app_jac_w_type, phi_op_type, A_type>;

  using base_pol_t = ::rompp::ode::policy::JacobianPolicyBase<this_t>;
  using base_incr_sol_t = rompp::rom::exp::IncrementalSolutionBase<this_t, app_state_w_type>;
  using scalar_type = typename core::details::traits<app_state_w_type>::scalar_t;
  
 private:
  mutable app_jac_w_type appJJ_;
  phi_op_type * phi_ = nullptr;
  A_type * A_ = nullptr;
  
  using base_incr_sol_t::y0FOM_;
  using base_incr_sol_t::yFOM_;
  
public:
  template <typename T=A_type,
   core::meta::enable_if_t<std::is_void<T>::value> * = nullptr>
  RomLSPGJacobianPolicy(const app_state_w_type & y0fom,
  			const app_jac_w_type & j0fom,
  			phi_op_type & phiOp)
    : base_incr_sol_t(y0fom), appJJ_(j0fom), phi_(&phiOp){}
  
  RomLSPGJacobianPolicy() = delete;
  ~RomLSPGJacobianPolicy() = default;
  //----------------------------------------------------------------
  
  template <typename ode_state_t,
	    typename ode_jac_t,
	    typename app_t>
  void operator()(const ode_state_t & odeY,
  		  ode_jac_t & odeJJ,
  		  const app_t & app,
  		  scalar_type t,
  		  scalar_type dt) const{
    // here y_n = odeY is reduced state
    // need to compute: Jac of R( phi y_n) = phi y_n - phi y_n-1 - dt * f(phi y)
    // let q be full state, then we have dR/dy = dR/dq dq/dy
    // here, we have: dR/dq stored into odeJJ
    
    // odeY is the REDUCED state, we need to reconstruct FOM state
    phi_->apply(odeY, yFOM_);
    
    // // since we are advancing the Incremental Solution,
    // // to compute the app residual we need to add the
    // // FOM initial condition to get full state 
    // yFOM_ += (*y0FOM_);
    
    // /// query the application for the jacobian
    // app.jacobian(*yFOM_.data(), *appJJ_.data(), t);

    // // do time discrete jacobian, which yields:
    // // appJJ = dR/dq 
    // ode::impl::implicit_euler_time_discrete_jacobian(appJJ_, dt);

    // // since dq/dy = phi_op
    // // right multiply appJJ with phi_op: odeJJ x phi
    // // this is = dR/dq dq/dy 
    // odeJJ = phi_->applyRight(appJJ_);
    
    // /// apply weighting
    // if (A_)
    //   A_->applyTranspose(appRHS_, odeR);
  }
  //----------------------------------------------------------------
  
  
  template <typename ode_state_t,
	    typename app_t>
  auto operator()(const ode_state_t & odeY,
  		  const app_t & app,
  		  scalar_type t,
  		  scalar_type dt) const{
    //   // odeY is the REDUCED state, we need to reconstruct FOM state
    phi_->apply(odeY, yFOM_);
    //   // since we are advancing the Incremental Solution,
    //   // to compute the app residual we need to add the
    //   // FOM initial condition to get full state 
    //   yFOM_ += (*y0FOM_);
    
    /// query the application for jacobian
    auto JJ = app.jacobian(*yFOM_.data(), t);
    core::Matrix<decltype(JJ)> JJw(JJ);
    
    //   // do time discrete residual
    //   ode::impl::implicit_euler_time_discrete_jacobian(JJw, dt);

    //   app_jac_w_type JJout(phi_->applyRight(JJw));
    return JJw;
    
    //   // /// apply weighting
    //   // if (A_)
    //   //   A_->applyTranspose(appRHS_, odeR);
  }
  //----------------------------------------------------------------


private:
  friend base_pol_t;
  friend base_incr_sol_t;

};//end class
  
}}}//end namespace rompp::rom::exp
#endif 
