
#ifndef ODE_IMPLICIT_POLICIES_META_HPP_
#define ODE_IMPLICIT_POLICIES_META_HPP_

#include "../base/ode_residual_policy_base.hpp"
#include "../base/ode_jacobian_policy_base.hpp"

#include "./euler/ode_implicit_euler_residual_standard_policy.hpp"
#include "./euler/ode_implicit_euler_jacobian_standard_policy.hpp"

namespace ode{
namespace meta {

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//-----------------------------------------------------------
// IMPLICIT EULER 
//-----------------------------------------------------------
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

template<typename policy_t, typename enable = void>
struct isLegitimateImplicitEulerResidualPolicy : std::false_type{};
  
template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type,
	  typename sizer_type,
	  typename ...Args>
struct isLegitimateImplicitEulerResidualPolicy<
  policy_t<state_type,residual_type,model_type,
  time_type,sizer_type,Args...>,
  typename std::enable_if<
    core::meta::publiclyInheritsFrom<
      policy_t<state_type,residual_type,model_type,
      time_type,sizer_type,Args...>,
      ode::policy::residualPolicyBase<policy_t,state_type,
						   residual_type,model_type,
						   time_type, sizer_type, Args...
						   >
                                    >::value
                         >::type
  > : std::true_type{};
//-----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isLegitimateImplicitEulerJacobianPolicy : std::false_type{};
  
template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type,
	  typename sizer_type,
	  typename ...Args>
struct isLegitimateImplicitEulerJacobianPolicy<
  policy_t<state_type,residual_type,model_type,
  time_type,sizer_type,Args...>,
  typename std::enable_if<
    core::meta::publiclyInheritsFrom<
      policy_t<state_type,residual_type,model_type,
      time_type,sizer_type,Args...>,
      ode::policy::jacobianPolicyBase<policy_t,state_type,
						   residual_type,model_type,
						   time_type, sizer_type, Args...
						   >
                                    >::value
                         >::type
  > : std::true_type{};
//-----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isImplicitEulerResidualStandardPolicy
  : std::false_type{};

template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type,
	  typename sizer_type>
struct isImplicitEulerResidualStandardPolicy<
  policy_t<state_type,residual_type,model_type,time_type, sizer_type>,
  typename std::enable_if<
    std::is_same<policy_t<state_type,residual_type,model_type,time_type, sizer_type>,
		 ode::policy::implicitEulerStandardResidual<
		   state_type, residual_type,model_type, time_type, sizer_type>
		 >::value
			  >::type
  > : std::true_type{};
//----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isImplicitEulerJacobianStandardPolicy
  : std::false_type{};

template <template <typename...> class policy_t,
	  typename state_type,
	  typename jacobian_type,
	  typename model_type,
	  typename time_type,
	  typename sizer_type>
struct isImplicitEulerJacobianStandardPolicy<
  policy_t<state_type,jacobian_type,model_type,time_type, sizer_type>,
  typename std::enable_if<
    std::is_same<policy_t<state_type,jacobian_type,model_type,time_type, sizer_type>,
		 ode::policy::implicitEulerStandardJacobian<
		   state_type, jacobian_type,model_type, time_type, sizer_type>
		 >::value
			  >::type
  > : std::true_type{};
//----------------------------------------------------------------


  
} // namespace meta
} // namespace core
#endif
