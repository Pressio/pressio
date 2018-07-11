
#ifndef ROM_GALERKIN_EXPLICIT_POLICY_HPP_
#define ROM_GALERKIN_EXPLICIT_POLICY_HPP_

#include "rom_ConfigDefs.hpp"
#include "ODE_POLICIES"

namespace rom{

// template<typename state_type,
// 	 typename residual_type,
// 	 typename model_type,
// 	 typename time_type,
// 	 typename lsv_type>
// class romGalerkinExplicit
//   : public ode::explicitResidualPolicyBase<
//   romGalerkinExplicit,
//   state_type, residual_type,
//   model_type, time_type, lsv_type>
// {
// public:
//   romGalerkinExplicit() = default;
//   ~romGalerkinExplicit() = default;

// private:
//   template <typename U = state_type,
// 	    typename T = residual_type,
// 	    typename
// 	    std::enable_if<
// 	      core::meta::is_coreVectorWrapper<U>::value==true &&
// 	      core::meta::is_coreVectorWrapper<T>::value==true
// 	      >::type * = nullptr
// 	    >
//   void computeImpl(const U & y, T & R,
// 		   model_type & model, time_type t){
//     model.residual(*y.data(), *R.data(), t);
//   }

// private:
//   friend explicitResidualPolicyBase<romGalerkinExplicit,
// 				    state_type, residual_type,
// 				    model_type, time_type, lsv_type
// 				    >;
// };//end class  

}//end namespace ode  
#endif 
