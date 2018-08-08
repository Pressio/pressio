
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
	 typename space_residual_type,
	 typename model_type,
	 typename phi_type,
	 typename A_type>
class RomGalerkinExplicitResidualPolicy
  : public ode::policy::ExplicitResidualPolicyBase<
  RomGalerkinExplicitResidualPolicy<state_type, space_residual_type,
				    model_type, phi_type, A_type>>
{

private:
  using this_t = RomGalerkinExplicitResidualPolicy<state_type,
						   space_residual_type,
						   model_type,
						   phi_type, A_type>;
  using base_t = ode::policy::ExplicitResidualPolicyBase<this_t>;
  
private:
  state_type yFOM_;
  space_residual_type appRHS_;
  phi_type * phi_;
  A_type * A_;
  
public:
  RomGalerkinExplicitResidualPolicy(const state_type & y0fom,
				    const space_residual_type & r0fom,
				    phi_type & phiOp,
				    A_type & AOp)
    : yFOM_(y0fom), appRHS_(r0fom), phi_(&phiOp), A_(&AOp)
  {}
  
  RomGalerkinExplicitResidualPolicy() = delete;
  ~RomGalerkinExplicitResidualPolicy() = default;  

private:
  template <typename T1 = state_type,
	    typename T2 = space_residual_type,
	    typename T3 = model_type,
	    typename scalar_type>
  void computeImpl(const T1 & y,
		   T2 & odeR,
		   T3 & model,
		   scalar_type t)
  {
    // y coming in is the REDUCED state, so we need to reconstruct full state
    assert( y.globalSize() == phi_->globalCols() );
    yFOM_ = core::matrixVectorProduct(*phi_, y);
    std::cout << "yFOM_ = " << yFOM_.globalSize() << std::endl;

    // compute space residual from the application
    model.residual(*yFOM_.data(), *appRHS_.data(), t);

    // apply weighting 
    odeR = A_->apply(appRHS_);
    std::cout << "odeR = " << odeR.globalSize() << std::endl;
  }
  
private:
  friend base_t;
  
};//end class
  
}//end namespace exp
}//end namespace rom
#endif 
