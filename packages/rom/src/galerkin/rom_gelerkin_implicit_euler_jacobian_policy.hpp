
#ifndef ROM_GALERKIN_IMPLICIT_EULER_JACOBIAN_POLICY_HPP_
#define ROM_GALERKIN_IMPLICIT_EULER_JACOBIAN_POLICY_HPP_

#include "../rom_forward_declarations.hpp"
#include "../../../CORE_ALL"
#include "../../../ode/src/implicit/ode_jacobian_impl.hpp"
#include "../../../ode/src/implicit/policies/base/ode_jacobian_policy_base.hpp"
// #include "../rom_incremental_solution_base.hpp"

namespace rompp{ namespace rom{ 

// template<typename app_state_w_type,
// 	 typename app_jac_w_type,
// 	 typename phi_op_type,  
// 	 typename ode_state_w_type,
// 	 typename ode_jac_w_type,
// 	 typename A_type>
// class RomGalerkinImplicitJacobianPolicy<ode::ImplicitEnum::Euler,
// 					app_state_w_type, app_jac_w_type,
// 					phi_op_type, ode_state_w_type,
// 					ode_jac_w_type, A_type>
//   : // inherits from jacobian policy base
//     public ode::policy::JacobianPolicyBase<
// 		RomGalerkinImplicitJacobianPolicy<
// 		  ode::ImplicitEnum::Euler,
// 		  app_state_w_type, app_jac_w_type,
// 		  phi_op_type, ode_state_w_type,
// 		  ode_jac_w_type, A_type>>,
//     // inherits from incremental solution
//     private IncrementalSolutionBase<
//                 RomGalerkinImplicitJacobianPolicy<
// 		  ode::ImplicitEnum::Euler,
// 		  app_state_w_type, app_jac_w_type,
// 		  phi_op_type, ode_state_w_type,
// 		  ode_jac_w_type, A_type>, app_state_w_type>
// {
  
//   using this_t = RomGalerkinImplicitJacobianPolicy<
//     ::rompp::ode::ImplicitEnum::Euler,
//     app_state_w_type, app_jac_w_type, phi_op_type,
//     ode_state_w_type,  ode_jac_w_type, A_type>;

//   using base_pol_t = ::rompp::ode::policy::JacobianPolicyBase<this_t>;
//   using base_incr_sol_t = IncrementalSolutionBase<this_t, app_state_w_type>;
//   using scalar_type = typename core::details::traits<app_state_w_type>::scalar_t;
  
//  private:
//   mutable app_jac_w_type appJJ_;
//   phi_op_type * phi_ = nullptr;
//   A_type * A_ = nullptr;  
//   using base_incr_sol_t::y0FOM_;
//   using base_incr_sol_t::yFOM_;
//   using app_jac_nat_type = typename core::details::traits<app_jac_w_type>::wrapped_t;
  
// public:

//   // default constructor deleted
//   RomGalerkinImplicitJacobianPolicy() = delete;
//   // default destructor
//   ~RomGalerkinImplicitJacobianPolicy() = default;
  
//   template <typename T=A_type,
//    ::rompp::mpl::enable_if_t<std::is_same<T,phi_op_type>::value> * = nullptr>
//   RomGalerkinImplicitJacobianPolicy(const app_state_w_type & y0fom,
// 					 const app_jac_w_type & j0fom,
// 					 phi_op_type & phiOp)
//     : base_incr_sol_t(y0fom), appJJ_(j0fom), phi_(&phiOp), A_(&phiOp){}

  
  
// public:
  
//   template <typename app_t>
//   ode_jac_w_type operator()(const ode_state_w_type & odeY,
// 			    const app_t & app,
// 			    scalar_type t,
// 			    scalar_type dt) const{
    
//     // need to compute I - dt * A^T J phi
//     reconstructFOMState(odeY);
    
//     /// query the application for the jacobian, i.e. J
//     app.jacobian(*yFOM_.data(), *appJJ_.data(), t);

//     // J phi
//     auto JJphi = phi_->applyRight(appJJ_);
//     // A^T J phi
//     ode_jac_w_type ATJJphi = A_->applyTranspose(JJphi);
    
//     // // compute time discrete residual
//     // ode::impl::implicit_time_discrete_jacobian(ATJJphi, dt);    
//     return ATJJphi;
//   }
//   //----------------------------------------------------------------
  
//   template <typename app_t>
//   void operator()(const ode_state_w_type & odeY,
//   		  ode_jac_w_type & odeJJ,
//   		  const app_t & app,
//   		  scalar_type t,
//   		  scalar_type dt) const{

//     odeJJ = (*this)(odeY, app, t, dt);
//     // reconstructFOMState(odeY);    

//     // /// query the application for the jacobian
//     // app.jacobian(*yFOM_.data(), *appJJ_.data(), t);

//     // // J phi
//     // auto JJphi = phi_->applyRight(appJJ_);
//     // // A^T J phi
//     // A_->applyTranspose(JJphi, odeJJ);
    
//     // // do time discrete jacobian
//     // ode::impl::implicit_euler_time_discrete_jacobian(odeJJ, dt);
//   }
//   //----------------------------------------------------------------

//  private:
//   void reconstructFOMState(const ode_state_w_type & odeY)const {

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
