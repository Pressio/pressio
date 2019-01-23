
#ifndef ROM_GALERKIN_EXPLICIT_RESIDUAL_POLICY_HPP_
#define ROM_GALERKIN_EXPLICIT_RESIDUAL_POLICY_HPP_

#include "../rom_ConfigDefs.hpp"
#include "../../../CORE_ALL"
#include "../../../ode/src/explicit/policies/ode_explicit_residual_policy_base.hpp"
//#include "../rom_incremental_solution_base.hpp"
// #include "../rom_data_base.hpp"

namespace rompp{ namespace rom{

// template<typename state_type,
// 	 typename space_res_type,
// 	 typename phi_op_type,
// 	 typename A_type = phi_op_type>
// class RomGalerkinExplicitResidualPolicy
//   : public ode::policy::ExplicitResidualPolicyBase<
//                 RomGalerkinExplicitResidualPolicy<
// 		  state_type, space_res_type, phi_op_type, A_type>>,
//     protected RomStateData<state_type, phi_op_type, 0>,
//     protected RomRHSData<space_res_type, 1>{

//   using this_t = RomGalerkinExplicitResidualPolicy<state_type,
// 						   space_res_type,
// 						   phi_op_type,
// 						   A_type>;

//   using base_pol_t = ode::policy::ExplicitResidualPolicyBase<this_t>;

//   using base_state_data_t = ::rompp::rom::RomStateData<state_type, phi_op_type, 0>;

//   using base_rhs_data_t = ::rompp::rom::RomRHSData<space_res_type, 1>;


// private:
//   A_type * A_ = nullptr;

//   using base_state_data_t::y0FOM_;
//   using base_state_data_t::yFOM_;
//   using base_state_data_t::phi_;
//   using base_rhs_data_t::appRHS_;

// public:

//   template <typename T=A_type,
//    core::meta::enable_if_t<
//      std::is_same<T,phi_op_type>::value
//      > * = nullptr>
//   RomGalerkinExplicitResidualPolicy(const state_type & y0fom,
// 				    const space_res_type & r0fom,
// 				    phi_op_type & phiOp)
//     : base_state_data_t(y0fom, phiOp),
//       base_rhs_data_t(r0fom),
//       A_(&phiOp){}

//   RomGalerkinExplicitResidualPolicy() = delete;

//   ~RomGalerkinExplicitResidualPolicy() = default;

// public:

//   template <typename ode_state_t,
// 	    typename ode_res_t,
// 	    typename app_type,
// 	    typename scalar_type>
//   void operator()(const ode_state_t & odeY,
// 		  ode_res_t & odeR,
// 		  const app_type & app,
// 		  scalar_type t) const {

//     /*- odeY: has type = ode_state_t, NOT for sure = state_type
//      *  - odeR: has type = ode_res_t, NOT for sure = space_res_type

//      * Typically ode_types are different!
//      *
//      * Because the way we compute ROM stuff does not have to
//      * be the same as the types used by the application.
//      * It can be, but not necessarily.
//      */

//     // odeY is the REDUCED state, we need to reconstruct FOM state
//     phi_->apply(odeY, yFOM_);

//     // since we are advancing the Incremental Solution,
//     // to compute the app residual we need to add the
//     // FOM initial condition to get full state
//     yFOM_ += *y0FOM_;

//     /// query the application for the SPACE residual
//     app.residual(*yFOM_.data(), *appRHS_[0].data(), t);

//     /// apply weighting
//     A_->applyTranspose(appRHS_[0], odeR);
//   }

// private:
//   friend base_pol_t;

// };//end class

}}//end namespace rompp::rom
#endif




// template <typename T = A_type,
// 	    core::meta::enable_if_t<
// 	      !std::is_void<A_type>::value,
// 	      T> * = nullptr
// 	    >
// RomGalerkinExplicitResidualPolicy(const state_type & y0fom,
// 				    const space_res_type & r0fom,
// 				    phi_op_type & phiOp,
// 				    T & AOp)
//   : base_incr_sol_t(y0fom), appRHS_(r0fom), phi_(&phiOp), A_(&AOp){}

// RomGalerkinExplicitResidualPolicy(const state_type & y0fom,
// 				    const space_res_type & r0fom,
// 				    phi_op_type & phiOp)
//   : base_incr_sol_t(y0fom), appRHS_(r0fom), phi_(&phiOp){}

// RomGalerkinExplicitResidualPolicy(const state_type & y0fom,
// 				    const space_res_type & r0fom)
//   : base_incr_sol_t(y0fom), appRHS_(r0fom){}
