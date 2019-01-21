
#ifndef ROM_LSPG_JACOBIAN_POLICY_HPP_
#define ROM_LSPG_JACOBIAN_POLICY_HPP_

#include "../rom_forward_declarations.hpp"
#include "../../../CORE_ALL"
#include "../../../ode/src/implicit/ode_jacobian_impl.hpp"
#include "../../../ode/src/implicit/policies/base/ode_jacobian_policy_base.hpp"
#include "../rom_data_base.hpp"

namespace rompp{ namespace rom{

template< typename app_state_w_type,
	  typename jac_type,
	  typename phi_op_type,
	  int maxNstates,
	  typename A_type>
class RomLSPGJacobianPolicy
  : public ode::policy::JacobianPolicyBase<
		RomLSPGJacobianPolicy<app_state_w_type,
				      jac_type,
				      phi_op_type,
				      maxNstates,
				      A_type>>,
    protected RomStateData<app_state_w_type, phi_op_type,  maxNstates>
{

  using this_t 		= RomLSPGJacobianPolicy<app_state_w_type,
						jac_type,
						phi_op_type,
						maxNstates,
						A_type>;

  using base_pol_t 	= ::rompp::ode::policy::JacobianPolicyBase<this_t>;

  using base_state_data_t = ::rompp::rom::RomStateData<app_state_w_type, phi_op_type, maxNstates>;

  using scalar_type 	= typename core::details::traits<app_state_w_type>::scalar_t;

  using base_state_data_t::phi_;
  using base_state_data_t::yFOM_;

public:

  RomLSPGJacobianPolicy() = delete;

  ~RomLSPGJacobianPolicy() = default;

  template <typename T=A_type,
	    core::meta::enable_if_t<std::is_void<T>::value> * = nullptr>
    RomLSPGJacobianPolicy(const app_state_w_type & y0fom,
			  phi_op_type & phiOp)
    : base_state_data_t(y0fom, phiOp){}


 public:

  template <::rompp::ode::ImplicitEnum odeMethod,
	     typename ode_state_t,
	     typename app_t>
  jac_type operator()(const ode_state_t & odeY,
		      const app_t & app,
		      scalar_type t,
		      scalar_type dt) const
  {
    base_state_data_t::template reconstructCurrentFOMState(odeY);

    // compute the Jac phi product, where Jac is the spatial jacobian
    auto * basis = phi_->getOperator();
    if (!JJ_){
      auto res = app.applyJacobian(*yFOM_.data(), *basis->data(), t);
      JJ_ = std::make_shared<jac_type>(res);
    } else{
      app.applyJacobian(*yFOM_.data(), *basis->data(), *JJ_->data(), t);
    }

    ode::impl::implicit_time_discrete_jacobian<odeMethod>(*JJ_, dt, *basis);

    app.applyPreconditioner(*yFOM_.data(), *JJ_->data(), t);

    // need to apply final weighting if any
    // if (A_) A_->apply....
    return *JJ_;
  }


  template <::rompp::ode::ImplicitEnum odeMethod,
	     typename ode_state_t,
	     typename ode_jac_t,
	     typename app_t>
  void operator()(const ode_state_t & odeY,
  		  ode_jac_t & odeJJ,
  		  const app_t & app,
  		  scalar_type t,
  		  scalar_type dt) const
  {
    base_state_data_t::template reconstructCurrentFOMState(odeY);

    auto * basis = phi_->getOperator();
    app.applyJacobian(*yFOM_.data(), *basis->data(), *odeJJ.data(), t);
    ode::impl::implicit_time_discrete_jacobian<odeMethod>(odeJJ, dt, *basis);
    app.applyPreconditioner(*yFOM_.data(), *odeJJ.data(), t);
  }


private:
  friend base_pol_t;

  mutable std::shared_ptr<jac_type> JJ_ = nullptr;

};//End class

}}//end namespace rompp::rom
#endif
