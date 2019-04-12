
#ifndef ROM_GALERKIN_IMPLICIT_EULER_RESIDUAL_POLICY_HPP_
#define ROM_GALERKIN_IMPLICIT_EULER_RESIDUAL_POLICY_HPP_

#include "../rom_forward_declarations.hpp"
#include "../../../CORE_ALL"
#include "../../../ode/src/implicit/ode_residual_impl.hpp"
#include "../../../ode/src/implicit/policies/base/ode_implicit_residual_policy_base.hpp"
// #include "../rom_incremental_solution_base.hpp"

namespace rompp{ namespace rom{ 
    
// template<typename app_state_w_type,
// 	 typename app_res_w_type,
// 	 typename phi_op_type,
// 	 typename ode_state_w_type,
// 	 typename ode_res_w_type,
// 	 typename A_type>
// class RomGalerkinImplicitResidualPolicy<ode::ImplicitEnum::Euler,
// 					app_state_w_type, app_res_w_type,
// 					phi_op_type, ode_state_w_type,
// 					ode_res_w_type, A_type>
//   : // inherits from implicit residual policy base
//     public ode::policy::ImplicitResidualPolicyBase<
//                 RomGalerkinImplicitResidualPolicy<
// 		  ode::ImplicitEnum::Euler,
// 		  app_state_w_type, app_res_w_type,
// 		  phi_op_type, ode_state_w_type,
// 		  ode_res_w_type, A_type>>,
//     // inherits from incremental solution
//     private IncrementalSolutionBase<
//                 RomGalerkinImplicitResidualPolicy<
// 		  ::rompp::ode::ImplicitEnum::Euler,
// 		  app_state_w_type, app_res_w_type,
// 		  phi_op_type, ode_state_w_type,
// 		  ode_res_w_type, A_type>, app_state_w_type>
// {
  
//   using this_t = RomGalerkinImplicitResidualPolicy<
//     ::rompp::ode::ImplicitEnum::Euler,
//     app_state_w_type, app_res_w_type, phi_op_type,
//     ode_state_w_type, ode_res_w_type, A_type>;

//   using base_pol_t = ode::policy::ImplicitResidualPolicyBase<this_t>;
//   using base_incr_sol_t = IncrementalSolutionBase<this_t, app_state_w_type>;
//   using scalar_type = typename core::details::traits<app_state_w_type>::scalar_t;
  
//  private:
//   mutable app_res_w_type appRHS_;
//   mutable app_state_w_type yFOMnm1_;
//   phi_op_type * phi_ = nullptr;
//   A_type * A_ = nullptr;  
//   using base_incr_sol_t::y0FOM_;
//   using base_incr_sol_t::yFOM_;
 

// public:

//   // default constructor deleted
//   RomGalerkinImplicitResidualPolicy() = delete;
//   // default destructor
//   ~RomGalerkinImplicitResidualPolicy() = default;

//   template <typename T=A_type,
//    ::rompp::mpl::enable_if_t<std::is_same<T,phi_op_type>::value> * = nullptr>
//   RomGalerkinImplicitResidualPolicy(const app_state_w_type & y0fom,
// 				    const app_res_w_type & r0fom,
// 				    phi_op_type & phiOp)
//     : base_incr_sol_t(y0fom), appRHS_(r0fom),
//       yFOMnm1_(y0fom), phi_(&phiOp), A_(&phiOp){}

  
// public:
  
//   //---------------------------------------------------------------------------
//   // compute: R( hat(y)_n) = hat(y)_n - hat(y)_n-1 - dt * A^T * f(phi hat(y))  
//   //---------------------------------------------------------------------------
//   template <typename app_t>
//   ode_res_w_type operator()(const ode_state_w_type & odeY,
// 			    const std::array<ode_state_w_type, 1> & prevYs,
// 			    const app_t & app,
// 			    scalar_type t,
// 			    scalar_type dt) const{
    
//     ode_res_w_type modRHS(odeY);//this is just to construct it 
//     modRHS.setZero();
//     // (*this)(odeY, modRHS, prevYs, app, t, dt);


//    // // odeY is the REDUCED state, we need to reconstruct FOM state
//    // reconstructFOMState(odeY, prevYs[0]);
//    // /// query the application for the SPACE residual 
//    // app.residual(*yFOM_.data(), *appRHS_.data(), t);
//    // /// apply operator
//    // A_->applyTranspose(appRHS_, modRHS);    
//    // /// do time discrete residual
//    // ode::impl::implicit_euler_time_discrete_residual(odeY, prevYs[0], modRHS, dt);
//     return modRHS;
//   }  
//   //----------------------------------------------------------------
  
//   template <typename app_t>
//   void operator()(const ode_state_w_type & odeY,
//   		  ode_res_w_type & odeR,
//   		  const std::array<ode_state_w_type, 1> & prevYs,
//   		  const app_t & app,
//   		  scalar_type t,
//   		  scalar_type dt) const{

//     // // odeY is the REDUCED state, we need to reconstruct FOM state
//     // reconstructFOMState(odeY, prevYs[0]);

//     // /// query the application for the SPACE residual 
//     // app.residual(*yFOM_.data(), *appRHS_.data(), t);

//     // /// apply operator
//     // A_->applyTranspose(appRHS_, odeR);

//     // /// do time discrete residual
//     // ode::impl::implicit_euler_time_discrete_residual(odeY, prevYs[0], odeR, dt);
//   }
//   //----------------------------------------------------------------

//  private:
//   void reconstructFOMState(const ode_state_w_type & odeY,
// 			   const ode_state_w_type & odeYm1)const {
//     // odeY is the REDUCED state, reconstruct FOM state and add FOM(t=0)
//     // because we are advancing the Incremental Solution
//     phi_->apply(odeY, yFOM_);
//     yFOM_ += (*y0FOM_);
//   }
//   //----------------------------------------------------------------
  
// private:
//   friend base_pol_t;
//   friend base_incr_sol_t;

// };//end class
  
}}//end namespace rompp::rom
#endif 
    
