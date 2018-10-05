
#ifndef ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HELPER_USING_ARB_POLICY_HPP_
#define ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HELPER_USING_ARB_POLICY_HPP_

#include "./impl/ode_explicit_euler_stepper_impl.hpp"
#include "ode_explicit_stepper_helper_info.hpp"

namespace rompp{ namespace ode{ namespace impl{
      
// ///////////////////////////////////////////////////
// /* Specialization for when we pass
//     state_type 
//     arbitrary_type
// Here we have: state_type = residual_type
// The arbitrary residual policy object itself
// needs to have inside the app object to compute */
// ///////////////////////////////////////////////////

// template<ExplicitSteppersEnum whichone,
// 	 typename T1, typename T2>
// class ExplicitStepperHelper<
//   whichone,
//   T1,T2,
//   core::meta::enable_if_t<
//     ode::meta::is_expl_type_pair_state_and_respol<T1,T2>::value
//     >
//   > : public
// stepper_helper_info<whichone,
//   typename ode::meta::is_expl_type_pair_state_and_respol<T1,T2>::state_type,
//   void,
//   typename ode::meta::is_expl_type_pair_state_and_respol<T1,T2>::state_type,
//   typename ode::meta::is_expl_type_pair_state_and_respol<T1,T2>::residual_policy_type
//   >::base_std_t{

//   using types_pair = ode::meta::is_expl_type_pair_state_and_respol<T1, T2>;
//   using ode_state_t = typename types_pair::state_type;
//   using res_pol_t = typename types_pair::residual_policy_type;
  
//   using info_t = stepper_helper_info<whichone, ode_state_t, void,
// 				     ode_state_t, res_pol_t>;
//   using pol_t = typename info_t::res_pol_t;
  
// public:  
//   using base_t = typename info_t::base_std_t;

// public:

//   template < typename... Args>
//   ExplicitStepperHelper(pol_t & policyObj,
// 			ode_state_t const & y0,
// 			ode_state_t const & r0,
// 			Args&&... rest)
//     : base_t(policyObj, y0, r0,
// 	     std::forward<Args>(rest)...){}

//   ExplicitStepperHelper() = delete;
//   ~ExplicitStepperHelper() = default;

// };//end class



// ///////////////////////////////////////////////////
// // Specialization for when we pass 
// //    state_type
// //    arbitrary_policy_type
// //    model_type 
// // Here we have: state_type = residual_type
// ///////////////////////////////////////////////////

// template<ExplicitSteppersEnum whichone,
// 	 typename T1, typename T2, typename T3>
// class ExplicitStepperHelper<
//   whichone,
//   T1,T2,T3,
//   core::meta::enable_if_t<
//     ode::meta::is_expl_type_set_state_respol_model<T1,T2,T3>::value
//     >
//   > : public
// stepper_helper_info<whichone,
//    typename ode::meta::is_expl_type_set_state_respol_model<T1,T2,T3>::state_type,
//    typename ode::meta::is_expl_type_set_state_respol_model<T1,T2,T3>::model_type,
//    typename ode::meta::is_expl_type_set_state_respol_model<T1,T2,T3>::state_type,
//    typename ode::meta::is_expl_type_set_state_respol_model<T1,T2,T3>::residual_policy_type
//   >::base_std_t{

//   using types_seq = ode::meta::is_expl_type_set_state_respol_model<T1,T2,T3>;
//   using ode_state_t = typename types_seq::state_type;
//   using model_t = typename types_seq::model_type;
//   using res_pol_t = typename types_seq::residual_policy_type;

//   using info_t = stepper_helper_info<whichone, ode_state_t, model_t,
// 				     ode_state_t, res_pol_t>;
//   using pol_t = typename info_t::res_pol_t;
  
// public:  
//   using base_t = typename info_t::base_std_t;

// public:

//   template < typename... Args>
//   ExplicitStepperHelper(model_t & modelObj,
// 			pol_t & policyObj,
// 			ode_state_t const & y0,
// 			ode_state_t const & r0,
// 			Args&&... rest)
//     : base_t(modelObj, policyObj, y0, r0,
// 	     std::forward<Args>(rest)...){}

//   ExplicitStepperHelper() = delete;
//   ~ExplicitStepperHelper() = default;

// };//end class

      
}}} // end namespace rompp::ode::impl
#endif 
