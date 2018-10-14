
#ifndef ROM_LSPG_RESIDUAL_POLICY_HPP_
#define ROM_LSPG_RESIDUAL_POLICY_HPP_

#include "rom_forward_declarations.hpp"
#include "../../../CORE_ALL"
#include "../../../ode/src/implicit/policies/base/ode_implicit_residual_policy_base.hpp"
#include "../../../ode/src/implicit/policies/base/ode_jacobian_policy_base.hpp"
#include "rom_incremental_solution_base.hpp"

namespace rompp{ namespace rom{ namespace exp{

// template<typename app_state_w_type,
// 	 typename app_res_w_type,
// 	 // basis operator type
// 	 typename phi_op_type,  
// 	 // weighting matrix by default is = void
// 	 typename A_type>
// class RomLSPGResidualPolicy<::rompp::ode::ImplicitSteppersEnum::Euler,
// 			     app_state_w_type, app_res_w_type,
// 			     phi_op_type, A_type>
//   : public ::rompp::ode::policy::ImplicitResidualPolicyBase<
//       RomLSPGResidualPolicy<::rompp::ode::ImplicitSteppersEnum::Euler,
// 			     app_state_w_type, app_res_w_type,
// 			     phi_op_type, A_type>, 1, 0>,
//     private IncrementalSolutionBase<
//       RomLSPGResidualPolicy<::rompp::ode::ImplicitSteppersEnum::Euler,
// 			     app_state_w_type, app_res_w_type,
// 			     phi_op_type, A_type>>{
  
//   using this_t = RomLSPGResidualPolicy<::rompp::ode::ImplicitSteppersEnum::Euler,
//     app_state_w_type, app_res_w_type, phi_op_type, A_type>;

//   using base_pol_t = ::rompp::ode::policy::ImplicitResidualPolicyBase<this_t>, 1, 0>;
//   using base_incr_sol_t = rompp::rom::exp::IncrementalSolutionBase<this_t, app_state_w_type>;
//   using scalar_type = typename core::details::traits<app_state_w_type>::scalar_t;
  
//  private:
//   mutable app_res_w_type appRHS_;
//   phi_op_type * phi_ = nullptr;
//   A_type * A_ = nullptr;
  
//   using base_incr_sol_t::y0FOM_;
//   using base_incr_sol_t::yFOM_;
  
// public:
//   template <typename T=A_type,
//    core::meta::enable_if_t<std::is_void<T>::value> * = nullptr>
//   RomLSPGResidualPolicy(const state_type & y0fom,
// 			const space_res_type & r0fom,
// 			phi_op_type & phiOp)
//     : base_incr_sol_t(y0fom), appRHS_(r0fom),
//       phi_(&phiOp), yFOMnm1_(y0fom){}
  
//   RomLSPGResidualPolicy() = delete;
//   ~RomLSPGResidualPolicy() = default;
//   //----------------------------------------------------------------
  
//   void operator()(const ode_state_type & odeY,
// 		  ode_res_type & odeR,
// 		  const std::array<ode_state_type, 1> & oldYs,
// 		  const model_type & app,
// 		  scalar_type t,
// 		  scalar_type dt) const{

//     // here y_n = odeY is reduced state
//     // need to compute: R( phi y_n) = phi y_n - phi y_n-1 - dt * f(phi y)
    
//     // odeY is the REDUCED state, we need to reconstruct FOM state
//     phi_->apply(odeY, yFOM_);

//     // reconstruct FOM state at previous step n-1
//     phi_->apply(oldYs[0], yFOMnm1_);
    
//     // since we are advancing the Incremental Solution,
//     // to compute the app residual we need to add the
//     // FOM initial condition to get full state 
//     yFOM_ += (*y0FOM_);
//     yFOMnm1_ += (*y0FOM_);
    
//     // /// query the application for the SPACE residual 
//     // app.residual(*yFOM_.data(), *appRHS_.data(), t);

//     // // do time discrete residual
//     // ::rompp::ode::impl::implicit_euler_time_discrete_residual(yFOM_, y0FOMnm1_,
//     // 							      R, dt);
        
//     // /// apply weighting
//     // if (A_)
//     //   A_->applyTranspose(appRHS_, odeR);
//   }
//   //----------------------------------------------------------------

// private:
//   friend base_pol_t;
//   friend base_incr_sol_t;
//   mutable state_type yFOMnm1_;

// };//end class
  
}}}//end namespace rompp::rom::exp
#endif 






  // // maybe sfinae here when types do not coincide
  // // with those above
  // template <typename ode_state_t,
  // 	    typename ode_res_t,
  // 	    typename app_type,
  // 	    typename scalar_type>
  // void operator()(const ode_state_t & odeY,
  // 		  ode_res_t & odeR,
  // 		  const app_type & app,
  // 		  scalar_type t) const {

  //   /*- odeY: has type = ode_state_t, NOT for sure = state_type
  //    *  - odeR: has type = ode_res_t, NOT for sure = space_res_type
  //    * 
  //    * Typically ode_types are different! 
  //    *
  //    * Because the way we compute ROM stuff does not have to 
  //    * be the same as the types used by the application. 
  //    * It can be, but not necessarily.
  //    */
    
  //   // odeY is the REDUCED state, we need to reconstruct FOM state
  //   phi_->apply(odeY, yFOM_);

  //   // since we are advancing the Incremental Solution,
  //   // to compute the app residual we need to add the
  //   // FOM initial condition to get full state 
  //   yFOM_ += *y0FOM_;
    
  //   /// query the application for the SPACE residual 
  //   app.residual(*yFOM_.data(), *appRHS_.data(), t);
    
  //   /// apply weighting
  //   A_->applyTranspose(appRHS_, odeR);
  // }
