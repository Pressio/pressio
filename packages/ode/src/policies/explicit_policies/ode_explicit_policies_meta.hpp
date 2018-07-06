
#ifndef ODE_EXPLICIT_POLICIES_META_HPP_
#define ODE_EXPLICIT_POLICIES_META_HPP_

#include "./base/ode_explicit_policy_base.hpp"
#include "./euler/ode_explicit_euler_standard_policy.hpp"
// #include "ode_explicit_runge_kutta4_standard_policy.hpp"


namespace ode{
namespace meta {
  
template<typename policy_t, typename enable = void>
struct isLegitimateExplicitResidualPolicy : std::false_type{};
  
template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type,
	  typename ...Args>
struct isLegitimateExplicitResidualPolicy<
  policy_t<state_type, residual_type,model_type, time_type, Args...>,
  typename std::enable_if<
    core::meta::publiclyInheritsFrom<
      policy_t<state_type,residual_type,model_type,time_type,Args...>,
      ode::policy::explicitResidualPolicyBase<policy_t,state_type,
					      residual_type,model_type,
					      time_type, Args...
					      >
                                    >::value
                         >::type
  > : std::true_type{};
//-----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isExplicitEulerResidualStandardPolicy: std::false_type{};

template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type>
struct isExplicitEulerResidualStandardPolicy<
  policy_t<state_type, residual_type,model_type, time_type>,
  typename std::enable_if<
    std::is_same<policy_t<state_type, residual_type,model_type, time_type>,
		 ode::policy::explicitEulerStandardResidual<
		   state_type, residual_type,model_type, time_type>
		 >::value
			  >::type
  > : std::true_type{};
//----------------------------------------------------------------

// template<typename policy_t, typename enable = void>
// struct isExplicitRungeKutta4ResidualStandardPolicy : std::false_type{};

// template <template <typename...> class policy_t,
// 	  typename state_type,
// 	  typename residual_type,
// 	  typename model_type,
// 	  typename time_type>
// struct isExplicitRungeKutta4ResidualStandardPolicy<
//   policy_t<state_type, residual_type,model_type, time_type>,
//   typename std::enable_if<
//     std::is_same<policy_t<state_type, residual_type,model_type, time_type>,
// 		 ode::policy::explicitRungeKutta4StandardResidual<
// 		   state_type, residual_type,model_type, time_type>
// 		 >::value
// 			  >::type
//   > : std::true_type{};
// //-----------------------------------------------------------------
  

} // namespace meta
} // namespace core
#endif






// /*
// A policy for residual for an explicit stepper has to meet:
// * it should inherit from the explicitResidualPolicyBase

// * it should contain a void PRIVATE method:  void computeImpl()
// if(a) is met, then if this is not there, then compiler gives an error

// * state_type and residual_type should not be simply std::is_integral
// * state_type and residual_type should not be simply std::is_floating_point

// * time_type should be std::is_floating_point

// **** more if comes to mind
// */
  
// template<typename policy_t, typename enable = void>
// struct isLegitimateOdeExplicitResidualPolicy : std::false_type{};
  
// template <template <typename...> class policy_t,
// 	  typename state_type,
// 	  typename residual_type,
// 	  typename model_type,
// 	  typename time_type,
// 	  typename ...Args>
// struct isLegitimateOdeExplicitResidualPolicy<
//   policy_t<state_type, residual_type,model_type, time_type, Args...>,
//   typename std::enable_if<
//     // first check that inheritance property
//     core::meta::publiclyInheritsFrom<
//       policy_t<state_type,residual_type,model_type,time_type,Args...>,
//       ode::policy::explicitResidualPolicyBase<policy_t,state_type,
// 					      residual_type,model_type,
// 					      time_type, Args...
// 					      >
//                                     >::value
//     // // check that state_type is legitimate
//     // isLegitimateStateType<state_type>::value &&
//     // // check that residual_type is legitimate 
//     // isLegitimateResidualType<state_type>::value &&
//     // // check that time_type is floating point
//     // std::is_floating_point<time_type>::value
//                          >::type
//   >
// {
//   static_assert(isLegitimateStateType<state_type>::value, "difdjfdfj");
//   static_assert(std::is_floating_point<time_type>::value, "ddfdfdsfdsfdifdjfdfj");
//   static constexpr bool value = true;
// };
