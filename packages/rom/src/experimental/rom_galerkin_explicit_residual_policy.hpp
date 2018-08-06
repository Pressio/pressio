
#ifndef ROM_GALERKIN_EXPLICIT_RESIDUAL_POLICY_HPP_
#define ROM_GALERKIN_EXPLICIT_RESIDUAL_POLICY_HPP_

#include "rom_ConfigDefs.hpp"
//#include "rom_incremental_solution_base.hpp"
#include "policies/base/ode_explicit_residual_policy_base.hpp"
//#include "ode_residual_impl.hpp"

namespace rom{
namespace exp{

template<typename state_type,
	 typename space_residual_type,
	 typename model_type,
	 typename sizer_type,
	 typename phi_type,
	 typename A_type>
class romGalerkinExplicitResidualPolicy
  : public ode::policy::ExplicitResidualPolicyBase<
  romGalerkinExplicitResidualPolicy<state_type, space_residual_type,
				    model_type, sizer_type,
				    phi_type, A_type>>
{

private:
  using this_t = romGalerkinExplicitResidualPolicy<state_type,
						   space_residual_type,
						   model_type,
						   sizer_type,
						   phi_type, A_type>;
  using base_t = ode::policy::ExplicitResidualPolicyBase<this_t>;
  
private:
  state_type yFOM_;
  space_residual_type appRHS_;
  phi_type * phiT_;
  A_type * A_;
  
public:
  romGalerkinExplicitResidualPolicy(const state_type & y0fom,
				    const space_residual_type & r0fom,
				    phi_type & phiTOp,
				    A_type & AOp)
    : yFOM_(y0fom), appRHS_(r0fom), phiT_(&phiTOp), A_(&AOp)
  {}
  
  romGalerkinExplicitResidualPolicy() = delete;
  ~romGalerkinExplicitResidualPolicy() = default;  

private:
  template <typename T1 = state_type,
	    typename T2 = space_residual_type,
	    typename T3 = model_type,
	    typename scalar_type>
  void computeImpl(const T1 & y,
		   T2 & R,
		   T3 & model,
		   scalar_type t)
  {

  }
  
private:
  friend base_t;
  
};//end class
  
}//end namespace exp
}//end namespace rom
#endif 
