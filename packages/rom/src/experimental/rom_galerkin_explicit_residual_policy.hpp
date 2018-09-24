
#ifndef ROM_GALERKIN_EXPLICIT_RESIDUAL_POLICY_HPP_
#define ROM_GALERKIN_EXPLICIT_RESIDUAL_POLICY_HPP_

#include "../rom_ConfigDefs.hpp"
#include "../../../CORE_ALL"
#include "../../../ode/src/policies/base/ode_explicit_residual_policy_base.hpp"
#include "rom_incremental_solution_base.hpp"

namespace rompp{
namespace rom{
namespace exp{

template<typename state_type,
	 typename space_res_type,
	 typename phi_op,
	 typename A_type = void>
class RomGalerkinExplicitResidualPolicy
  : public ode::policy::ExplicitResidualPolicyBase<
                RomGalerkinExplicitResidualPolicy<
		  state_type, space_res_type, phi_op, A_type>>,
    private IncrementalSolutionBase<
	        RomGalerkinExplicitResidualPolicy<
                  state_type, space_res_type,phi_op, A_type>, state_type>{
  
  using this_t = RomGalerkinExplicitResidualPolicy<state_type,
						   space_res_type,
						   phi_op,
						   A_type>;
  using base_pol_t = ode::policy::ExplicitResidualPolicyBase<this_t>;
  using base_incr_sol_t = rompp::rom::exp::IncrementalSolutionBase<this_t, state_type>;
  
public:
  space_res_type appRHS_;
  phi_op * phi_;
  A_type * A_;

protected:
  using base_incr_sol_t::y0FOM_;
  using base_incr_sol_t::yFOM_;
  
public:
  template <typename T = A_type,
	    core::meta::enable_if_t<
	      !std::is_void<A_type>::value,
	      T> * = nullptr
	    >
  RomGalerkinExplicitResidualPolicy(const state_type & y0fom,
  				    const space_res_type & r0fom,
  				    phi_op & phiOp,
				    T & AOp)
    : base_incr_sol_t(y0fom), appRHS_(r0fom), phi_(&phiOp), A_(&AOp){}

  RomGalerkinExplicitResidualPolicy(const state_type & y0fom,
  				    const space_res_type & r0fom,
  				    phi_op & phiOp)
    : base_incr_sol_t(y0fom), appRHS_(r0fom), phi_(&phiOp){}

  RomGalerkinExplicitResidualPolicy(const state_type & y0fom,
				    const space_res_type & r0fom)
    : base_incr_sol_t(y0fom), appRHS_(r0fom){}
  
  RomGalerkinExplicitResidualPolicy() = delete;
  ~RomGalerkinExplicitResidualPolicy() = default;  

private:

  // maybe sfinae here when types do not coincide
  // with those above
  template <typename ode_state_t,
	    typename ode_res_t,
	    typename app_type,
	    typename scalar_type>
  void operator()(const ode_state_t & odeY,
		  ode_res_t & odeR,
		  app_type & app,
		  scalar_type t){

    // I need to fix a few things like not creating new data structures,
    // but the overall idea is here.
    
    /*- odeY: has type = ode_state_t, NOT for sure = state_type
      - odeR: has type = ode_res_t, NOT for sure = space_res_type

      Because the way we compute ROM stuff does not have to 
      be the same as the types used by the application. 
      It can be, but not necessarily.
     */
    
    // odeY is the REDUCED state, we need to reconstruct FOM state
    yFOM_ = phi_->apply(odeY);

    // we advance incremental solution: add initial condition to get state 
    yFOM_ += y0FOM_;
    
    /// compute space residual from the application
    app.residual(*yFOM_.data(), *appRHS_.data(), t);
    
    ///apply weighting
    auto res = A_->applyTransp(appRHS_);
    
    // obviously this needs to be fixed
    for (size_t i=0; i<res.size(); i++)
      odeR[i] = res[i];
  }

private:
  friend base_pol_t;
  friend base_incr_sol_t;
  
};//end class
  
}//end namespace exp
}//end namespace rom
}//end namespace rompp
#endif 
