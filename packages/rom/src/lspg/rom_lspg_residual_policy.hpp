
#ifndef ROM_LSPG_RESIDUAL_POLICY_HPP_
#define ROM_LSPG_RESIDUAL_POLICY_HPP_

#include "../rom_forward_declarations.hpp"
#include "../../../CORE_ALL"
#include "../../../ode/src/implicit/ode_residual_impl.hpp"
#include "../../../ode/src/implicit/policies/base/ode_implicit_residual_policy_base.hpp"
#include "../rom_incremental_solution_base.hpp"

namespace rompp{ namespace rom{

template< typename app_state_w_type,
	  typename app_res_w_type,
	  typename phi_op_type,
	  typename A_type>
class RomLSPGResidualPolicy
  : public ode::policy::ImplicitResidualPolicyBase<
                RomLSPGResidualPolicy<app_state_w_type, app_res_w_type,
				      phi_op_type, A_type>>,
    private IncrementalSolutionBase<
	       RomLSPGResidualPolicy<app_state_w_type, app_res_w_type,
				     phi_op_type, A_type>, app_state_w_type>
{

  using this_t 		= RomLSPGResidualPolicy<app_state_w_type, app_res_w_type, phi_op_type, A_type>;
  using base_pol_t 	= ::rompp::ode::policy::ImplicitResidualPolicyBase<this_t>;
  using base_incr_sol_t = ::rompp::rom::IncrementalSolutionBase<this_t, app_state_w_type>;
  using scalar_type 	= typename core::details::traits<app_state_w_type>::scalar_t;

  static constexpr std::size_t storageSize_ = 4;

 public:

  template <typename T=A_type,
	    core::meta::enable_if_t<std::is_void<T>::value> * = nullptr>
  RomLSPGResidualPolicy(const app_state_w_type & y0fom,
			const app_res_w_type & r0fom,
			phi_op_type & phiOp)
    : base_incr_sol_t(y0fom),
      appRHS_{r0fom, r0fom, r0fom, r0fom},
      phi_(&phiOp),
      yFOMold_{y0fom, y0fom, y0fom, y0fom}{}

  template <typename T=A_type,
	    core::meta::enable_if_t<std::is_void<T>::value> * = nullptr>
  RomLSPGResidualPolicy(const app_state_w_type & y0fom,
  			phi_op_type & phiOp)
    : base_incr_sol_t(y0fom),
      phi_(&phiOp),
      yFOMold_{y0fom, y0fom, y0fom, y0fom}{}

  RomLSPGResidualPolicy() = delete;
  ~RomLSPGResidualPolicy() = default;


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
    // odeY, oldYsa are REDUCED states, reconstruct FOM states
    reconstructFOMStates<numAuxStates>(odeY, oldYs);

    /// query the application for the SPACE residual
    app.residual(*yFOM_.data(), *appRHS_[0].data(), t);

    /// do time discrete residual
    ode::impl::implicit_time_discrete_residual<odeMethod,
    					       storageSize_>(yFOM_,
							     yFOMold_,
							     appRHS_[0],
							     dt);

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
    reconstructFOMStates<numAuxStates>(odeY, oldYs);

    app.residual(*yFOM_.data(), *odeR.data(), t);

    ode::impl::implicit_time_discrete_residual<odeMethod,
    					       storageSize_>(yFOM_,
							     yFOMold_,
							     odeR,
							     dt);

    // /// apply weighting
    // if (A_) A_->applyTranspose(appRHS_, odeR);
  }
  //----------------------------------------------------------------

 private:

  template <int numAuxStates,
	    typename ode_state_t,
            core::meta::enable_if_t<numAuxStates==1> * = nullptr>
  void reconstructFOMStates(const ode_state_t & odeY,
          const std::array<ode_state_t,numAuxStates> & odeYprev) const{

    phi_->apply(odeY, yFOM_);
    phi_->apply(odeYprev[0], yFOMold_[0]);
    yFOM_ += (*y0FOM_);
    yFOMold_[0] += (*y0FOM_);
  }
  //----------------------------------------------------------------

  template <int numAuxStates,
	    typename ode_state_t,
            core::meta::enable_if_t<numAuxStates==2> * = nullptr>
  void reconstructFOMStates(const ode_state_t & odeY,
          const std::array<ode_state_t,numAuxStates> & odeYprev) const{

    phi_->apply(odeY, yFOM_);
    phi_->apply(odeYprev[0], yFOMold_[0]);
    phi_->apply(odeYprev[1], yFOMold_[1]);
    yFOM_ += (*y0FOM_);
    yFOMold_[0] += (*y0FOM_);
    yFOMold_[1] += (*y0FOM_);
  }
  //----------------------------------------------------------------

 private:
  friend base_pol_t;
  friend base_incr_sol_t;

  using base_incr_sol_t::y0FOM_;
  using base_incr_sol_t::yFOM_;

  mutable std::array<app_res_w_type, storageSize_> appRHS_     = {};
  mutable std::array<app_state_w_type, storageSize_> yFOMold_  = {};
  phi_op_type * phi_                                = nullptr;
  A_type * A_                                       = nullptr;


};//end class

}}//end namespace rompp::rom
#endif
