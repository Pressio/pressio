
#ifndef ROM_LSPG_BDF2_JACOBIAN_POLICY_HPP_
#define ROM_LSPG_BDF2_JACOBIAN_POLICY_HPP_

#include "../rom_forward_declarations.hpp"
#include "../../../CORE_ALL"
#include "../../../ode/src/implicit/ode_jacobian_impl.hpp"
#include "../../../ode/src/implicit/policies/base/ode_jacobian_policy_base.hpp"
#include "../rom_incremental_solution_base.hpp"

namespace rompp{ namespace rom{

template<typename app_state_w_type,
	 typename jac_type,
	 typename phi_op_type,
	 typename A_type /*defailt = void*/>
class RomLSPGJacobianPolicy<::rompp::ode::ImplicitSteppersEnum::BDF2,
			     app_state_w_type, jac_type,
			     phi_op_type, A_type>
  : public ode::policy::JacobianPolicyBase<
                RomLSPGJacobianPolicy<::rompp::ode::ImplicitSteppersEnum::BDF2,
				       app_state_w_type, jac_type,
				       phi_op_type, A_type>>,
    private IncrementalSolutionBase<
                RomLSPGJacobianPolicy<::rompp::ode::ImplicitSteppersEnum::BDF2,
				       app_state_w_type, jac_type,
				       phi_op_type, A_type>, app_state_w_type>{

  using this_t 		= RomLSPGJacobianPolicy<::rompp::ode::ImplicitSteppersEnum::BDF2,
   			  app_state_w_type, jac_type, phi_op_type, A_type>;
  using base_pol_t 	= ::rompp::ode::policy::JacobianPolicyBase<this_t>;
  using base_incr_sol_t = rompp::rom::IncrementalSolutionBase<this_t, app_state_w_type>;
  using scalar_type 	= typename core::details::traits<app_state_w_type>::scalar_t;

 private:
  mutable std::shared_ptr<jac_type> JJ_ = nullptr;
  phi_op_type * phi_ 			= nullptr;
  A_type * A_ 				= nullptr;

  using base_incr_sol_t::y0FOM_;
  using base_incr_sol_t::yFOM_;

public:
  template <typename T=A_type,
   core::meta::enable_if_t<std::is_void<T>::value> * = nullptr>
  RomLSPGJacobianPolicy(const app_state_w_type & y0fom,
  			phi_op_type & phiOp)
    : base_incr_sol_t(y0fom), phi_(&phiOp){}

  RomLSPGJacobianPolicy() = delete;
  ~RomLSPGJacobianPolicy() = default;
  //----------------------------------------------------------------

  // let q be full state, then we have dR/dy = dR/dq dq/dy

  template <typename ode_state_t,
  	    typename app_t>
  jac_type operator()(const ode_state_t & odeY,
		      const app_t & app,
		      scalar_type t,
		      scalar_type dt) const
  {
    reconstructFOMState(odeY);

    // compute the Jac phi product, where Jac is the spatial jacobian
    auto * basis = phi_->getOperator();
    if (!JJ_){
      auto res = app.applyJacobian(*yFOM_.data(), *basis->data(), t);
      JJ_ = std::make_shared<jac_type>(res);
    } else{
      app.applyJacobian(*yFOM_.data(), *basis->data(), *JJ_->data(), t);
    }

    ode::impl::implicit_bdf2_time_discrete_jacobian(*JJ_, dt, *basis);

    // need to apply final weighting if any
    // if (A_) A_->apply....
    return *JJ_;
  }
  //----------------------------------------------------------------

  template <typename ode_state_t,
  	    typename ode_jac_t,
  	    typename app_t>
  void operator()(const ode_state_t & odeY,
  		  ode_jac_t & odeJJ,
  		  const app_t & app,
  		  scalar_type t,
  		  scalar_type dt) const
  {
    reconstructFOMState(odeY);
    auto * basis = phi_->getOperator();
    app.applyJacobian(*yFOM_.data(), *basis->data(), *odeJJ.data(), t);
    ode::impl::implicit_bdf2_time_discrete_jacobian(odeJJ, dt, *basis);
  }
  //----------------------------------------------------------------

private:
  template <typename ode_state_t>
  void reconstructFOMState(const ode_state_t & odeY)const {
    phi_->apply(odeY, yFOM_);
    yFOM_ += (*y0FOM_);
  }
  //----------------------------------------------------------------

private:
  friend base_pol_t;
  friend base_incr_sol_t;

};//end class

}}//end namespace rompp::rom
#endif
