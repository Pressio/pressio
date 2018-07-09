
#ifndef ODE_POLICIES_META_HPP_
#define ODE_POLICIES_META_HPP_

#include "../base/ode_residual_policy_base.hpp"
#include "../base/ode_jacobian_policy_base.hpp"
#include "../standard/ode_residual_standard_policy.hpp"
#include "../standard/ode_jacobian_standard_policy.hpp"

namespace ode{
namespace meta {
  
template<typename policy_t,
	 typename enable = void>
struct isLegitimateResidualPolicy : std::false_type{};
  
template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type,
	  typename sizer_type,
	  typename ...Args>
struct isLegitimateResidualPolicy<
  policy_t<state_type, residual_type,
	   model_type, time_type, sizer_type, Args...>,
  typename std::enable_if<
    core::meta::publiclyInheritsFrom<
      policy_t<state_type, residual_type,
	       model_type, time_type, sizer_type, Args...>,
      ode::policy::residualPolicyBase<policy_t, state_type,
					      residual_type, model_type,
					      time_type, sizer_type, Args...
					      >
                                    >::value
                         >::type
  > : std::true_type{};
//-----------------------------------------------------------------

  
template<typename policy_t,
	 typename enable = void>
struct isResidualStandardPolicy: std::false_type{};

template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type,
	  typename sizer_type>
struct isResidualStandardPolicy<
  policy_t<state_type, residual_type,
	   model_type, time_type, sizer_type>,
  typename std::enable_if<
    std::is_same<policy_t<state_type, residual_type,
			  model_type, time_type, sizer_type>,
		 ode::policy::residualStandardPolicy<
		   state_type, residual_type,
		   model_type, time_type, sizer_type>
		 >::value
			  >::type
  > : std::true_type{};
//----------------------------------------------------------------


template<typename policy_t,
	 typename enable = void>
struct isLegitimateJacobianPolicy : std::false_type{};
  
template <template <typename...> class policy_t,
	  typename state_type,
	  typename jacobian_type,
	  typename model_type,
	  typename time_type,
	  typename sizer_type,
	  typename ...Args>
struct isLegitimateJacobianPolicy<
  policy_t<state_type, jacobian_type,
	   model_type, time_type, sizer_type, Args...>,
  typename std::enable_if<
    core::meta::publiclyInheritsFrom<
      policy_t<state_type, jacobian_type,
	       model_type, time_type, sizer_type, Args...>,
      ode::policy::jacobianPolicyBase<policy_t, state_type,
					      jacobian_type, model_type,
					      time_type, sizer_type, Args...
					      >
                                    >::value
                         >::type
  > : std::true_type{};
//-----------------------------------------------------------------

template<typename policy_t,
	 typename enable = void>
struct isJacobianStandardPolicy: std::false_type{};

template <template <typename...> class policy_t,
	  typename state_type,
	  typename jacobian_type,
	  typename model_type,
	  typename time_type,
	  typename sizer_type>
struct isJacobianStandardPolicy<
  policy_t<state_type, jacobian_type,
	   		model_type, time_type, sizer_type>,
  typename std::enable_if<
    std::is_same<policy_t<state_type, jacobian_type,
			  model_type, time_type, sizer_type>,
		 ode::policy::jacobianStandardPolicy<
		   state_type, jacobian_type,
		   model_type, time_type, sizer_type>
		 >::value
			  >::type
  > : std::true_type{};
//----------------------------------------------------------------

  
} // namespace meta
} // namespace core
#endif
