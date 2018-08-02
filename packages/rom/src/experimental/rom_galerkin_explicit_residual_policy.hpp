
#ifndef ROM_GALERKIN_EXPLICIT_RESIDUAL_POLICY_HPP_
#define ROM_GALERKIN_EXPLICIT_RESIDUAL_POLICY_HPP_

#include "rom_ConfigDefs.hpp"
#include "rom_incremental_solution_base.hpp"
#include "policies/base/ode_explicit_residual_policy_base.hpp"
//#include "ode_residual_impl.hpp"

namespace rom{
namespace exp{

template<typename state_type,
	 typename space_residual_type,
	 typename model_type,
	 typename sizer_type>
class romGalerkinExplicitResidualPolicy
  : public ode::policy::ExplicitResidualPolicyBase<
  romGalerkinExplicitResidualPolicy<state_type, space_residual_type,
				    model_type, sizer_type>>
{

private:
  using this_t = romGalerkinExplicitResidualPolicy<state_type, space_residual_type,
						   model_type, sizer_type>;
  using base_t = ode::policy::ExplicitResidualPolicyBase<this_t>;
  
private:
  // phi_op_t * phiOp_;
  // wei_op_t * WOp_;
  space_residual_type appRHS_;
  state_type yFOM_;
  
public:
  romGalerkinExplicitResidualPolicy(){}
  ~romGalerkinExplicitResidualPolicy() = default;  

// private:
//   template <typename U = state_type,
// 	    typename T = space_residual_type,
// 	    typename std::enable_if<
// 	      core::meta::is_core_vector<U>::value==true &&
// 	      core::meta::is_core_vector<T>::value==true
// 	    >::type * = nullptr>
//   void computeImpl(const U & y,
// 		   T & R,
// 		   const std::array<U, 1> & oldYs,
// 		   model_type & model,
// 		   time_type t,
// 		   time_type dt)
//   {
    
//   }  

private:
  friend base_t;
  
};//end class
  
}//end namespace exp
}//end namespace rom
#endif 





  // romGalerkinExplicitResidualPolicy(const state_type & y0fom,
  // 				    const state_type & y0r,
  // 				    phi_op_t & phiOp,
  // 				    wei_op_t & WOp)
  //   : phiOp_(&phiOp), WOp_(&WOp), appRHS_(y0fom.size()), yFOM_(y0fom.size())
  // {}
