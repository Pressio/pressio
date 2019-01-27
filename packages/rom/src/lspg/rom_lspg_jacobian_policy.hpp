
#ifndef ROM_LSPG_JACOBIAN_POLICY_HPP_
#define ROM_LSPG_JACOBIAN_POLICY_HPP_

#include "../rom_forward_declarations.hpp"
#include "../../../ode/src/implicit/ode_jacobian_impl.hpp"
#include "../../../ode/src/implicit/policies/base/ode_jacobian_policy_base.hpp"
#include "../rom_data_fom_states.hpp"

namespace rompp{ namespace rom{

template< typename fom_states_data,
	  typename apply_jac_return_type,
	  typename fom_apply_jac_policy>
class LSPGJacobianPolicy
  : public ode::policy::JacobianPolicyBase<
	LSPGJacobianPolicy<fom_states_data,
			   apply_jac_return_type,
			   fom_apply_jac_policy>>,
    protected fom_states_data,
    protected fom_apply_jac_policy{

protected:
  using this_t = LSPGJacobianPolicy<fom_states_data,
				    apply_jac_return_type,
				    fom_apply_jac_policy>;
  friend ode::policy::JacobianPolicyBase<this_t>;

  using fom_states_data::yFom_;
  using fom_states_data::yFomOld_;
  using fom_states_data::maxNstates_;
  using fom_states_data::decoderObj_;

public:
  static constexpr bool isResidualPolicy_ = false;
  using apply_jac_return_t = apply_jac_return_type;

public:
  LSPGJacobianPolicy() = delete;
  ~LSPGJacobianPolicy() = default;
  LSPGJacobianPolicy(const fom_states_data & fomStates,
		     const fom_apply_jac_policy & applyJacFunctor,
		     const apply_jac_return_type & applyJacObj)
    : fom_states_data(fomStates),
      fom_apply_jac_policy(applyJacFunctor),
      JJ_(applyJacObj){}

 public:

  template <ode::ImplicitEnum odeMethod,
	    typename ode_state_t,
	    typename ode_jac_t,
	    typename app_t,
	    typename scalar_t>
  void operator()(const ode_state_t & odeY,
		  ode_jac_t	    & odeJJ,
  		  const app_t	    & app,
		  scalar_t	    t,
		  scalar_t	    dt) const
  {
    fom_states_data::template reconstructCurrentFomState(odeY);
    const auto & basis = decoderObj_.getJacobianRef();
    fom_apply_jac_policy::evaluate(app, yFom_, basis, odeJJ, t);
    ode::impl::time_discrete_jacobian<odeMethod>(odeJJ, dt, basis);
  }

  template <ode::ImplicitEnum odeMethod,
	    typename ode_state_t,
	    typename app_t,
	    typename scalar_t>
  apply_jac_return_t operator()(const ode_state_t & odeY,
				const app_t	     & app,
				scalar_t	     t,
				scalar_t	     dt) const
  {
    (*this).template operator()<odeMethod>(odeY, JJ_, app, t, dt);
    return JJ_;
  }

private:
  mutable apply_jac_return_t JJ_ = {};
};

}}//end namespace rompp::rom
#endif





// // compute the Jac phi product, where Jac is the spatial jacobian
// auto * basis = phi_->getOperator();
// if (!JJ_){
//   auto res = fom_apply_jac_policy::evaluate(app, yFOM_, *basis, t);
//   JJ_ = std::make_shared<apply_jac_return_type>(res);
// } else{
//   fom_apply_jac_policy::evaluate(app, yFOM_, *basis, *JJ_, t);
//   //app.applyJacobian(*yFOM_.data(), *basis->data(), *JJ_->data(), t);
// }
