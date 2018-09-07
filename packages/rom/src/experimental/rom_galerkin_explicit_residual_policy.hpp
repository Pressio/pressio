
#ifndef ROM_GALERKIN_EXPLICIT_RESIDUAL_POLICY_HPP_
#define ROM_GALERKIN_EXPLICIT_RESIDUAL_POLICY_HPP_

#include "rom_ConfigDefs.hpp"
#include "CORE_ALL"
//#include "rom_incremental_solution_base.hpp"
#include "policies/base/ode_explicit_residual_policy_base.hpp"
//#include "ode_residual_impl.hpp"

namespace rom{
namespace exp{

template<typename state_type,
	 typename space_res_type,
	 typename phi_op,
	 typename A_type = void>
class RomGalerkinExplicitResidualPolicy
  : public ode::policy::ExplicitResidualPolicyBase<
  RomGalerkinExplicitResidualPolicy<state_type, space_res_type,
				    phi_op, A_type>>
{

  using this_t = RomGalerkinExplicitResidualPolicy<state_type,
						   space_res_type,
						   phi_op,
						   A_type>;
  using base_t = ode::policy::ExplicitResidualPolicyBase<this_t>;
  
  state_type yFOM_;
  space_res_type appRHS_;
  phi_op * phi_;
  A_type * A_;
  
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
    : yFOM_(y0fom), appRHS_(r0fom), phi_(&phiOp), A_(&AOp){}

  RomGalerkinExplicitResidualPolicy(const state_type & y0fom,
  				    const space_res_type & r0fom,
  				    phi_op & phiOp)
    : yFOM_(y0fom), appRHS_(r0fom), phi_(&phiOp){}

  RomGalerkinExplicitResidualPolicy(const state_type & y0fom,
				    const space_res_type & r0fom)
    : yFOM_(y0fom), appRHS_(r0fom){}
  
  RomGalerkinExplicitResidualPolicy() = delete;
  ~RomGalerkinExplicitResidualPolicy() = default;  

private:

  // maybe sfinae here when types do not coincide
  // with those above
  template <typename ode_state_t,
	    typename ode_res_t,
	    typename app_type,
	    typename scalar_type>
  void computeImpl(const ode_state_t & odeY,
		   ode_res_t & odeR,
		   app_type & app,
		   scalar_type t)
  {
    /* types: 
       - y: has type ode_state_t which is NOT for sure same as state_type
       - odeR: has type ode_res_t which is NOT for sure same as space_res_type

       This is because the way we compute ROM stuff does not have to 
       be the same as the types used by the application. 
       It can be, but it is not.
     */

    // odeY is the REDUCED state, we need to reconstruct full state
    phi_->apply(odeY, yFOM_);

    // // compute space residual from the application
    // app.residual(*yFOM_.data(), *appRHS_.data(), t);
    
    // // apply weighting 
    // A_->apply(appRHS_, odeR);

    //////////////////
    
    // // y coming in is the REDUCED state, so we need to reconstruct full state
    // assert( y.globalSize() == phi_->globalCols() );
    // yFOM_ = core::matrixVectorProduct(*phi_, y);
    // std::cout << "yFOM_ = " << yFOM_.globalSize() << std::endl;

    // // compute space residual from the application
    // model.residual(*yFOM_.data(), *appRHS_.data(), t);

    // // apply weighting 
    // odeR = A_->apply(appRHS_);
    // std::cout << "odeR = " << odeR.globalSize() << std::endl;
  }
  
private:
  friend base_t;
  
};//end class
  
}//end namespace exp
}//end namespace rom
#endif 
