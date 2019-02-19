
#ifndef ROM_DEFAULT_GALERKIN_EXPLICIT_RESIDUAL_POLICY_HPP_
#define ROM_DEFAULT_GALERKIN_EXPLICIT_RESIDUAL_POLICY_HPP_

#include "../rom_forward_declarations.hpp"
#include "../../../ode/src/explicit/policies/ode_explicit_residual_policy_base.hpp"
#include "../rom_data_fom_states.hpp"

namespace rompp{ namespace rom{

template <typename fom_states_data,
	  typename fom_rhs_data,
	  typename decoder_t>
class DefaultGalerkinExplicitResidualPolicy
  : public ode::policy::ExplicitResidualPolicyBase<
       DefaultGalerkinExplicitResidualPolicy<fom_states_data,
					     fom_rhs_data,
					     decoder_t>>,
    protected fom_states_data,
    protected fom_rhs_data{

protected:
  using this_t = DefaultGalerkinExplicitResidualPolicy<
		      fom_states_data, fom_rhs_data, decoder_t>;
  friend ode::policy::ImplicitResidualPolicyBase<this_t>;

  const decoder_t & decoder_;
  using fom_states_data::yFom_;
  using fom_rhs_data::fomRhs_;

public:
  static constexpr bool isResidualPolicy_ = true;

public:
  DefaultGalerkinExplicitResidualPolicy() = delete;
  ~DefaultGalerkinExplicitResidualPolicy() = default;
  DefaultGalerkinExplicitResidualPolicy(const fom_states_data & fomStates,
					const fom_rhs_data & fomResids,
					const decoder_t & decoder)
    : fom_states_data(fomStates),
      fom_rhs_data(fomResids),
      decoder_(decoder){}

public:
  template <typename galerkin_state_t,
	    typename galerkin_residual_t,
	    typename fom_t,
	    typename scalar_t>
  void operator()(const galerkin_state_t  & romY,
		  galerkin_residual_t	  & romR,
  		  const fom_t		  & app,
		  scalar_t		  t) const
  {
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("galerkin explicit residual");
#endif

    fom_states_data::template reconstructCurrentFomState(romY);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("fom eval rhs");
    app.residual(*yFom_.data(), *fomRhs_.data(), t);
    timer->stop("fom eval rhs");
#else
    app.residual(*yFom_.data(), *fomRhs_.data(), t);
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("phi^T*fom_rhs");
    core::ops::dot(decoder_.getJacobianRef(), fomRhs_, romR);
    timer->stop("phi^T*fom_rhs");
#else
    core::ops::dot(decoder_.getJacobianRef(), fomRhs_, romR);
#endif

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("galerkin explicit residual");
#endif
  }

};//end class










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
