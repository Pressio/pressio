
#ifndef ODE_IMPLICIT_POLICIES_META_HPP_
#define ODE_IMPLICIT_POLICIES_META_HPP_

#include "./base/ode_implicit_euler_residual_policy_base.hpp"
#include "./base/ode_implicit_euler_jacobian_policy_base.hpp"
#include "./euler/ode_implicit_euler_residual_standard_policy.hpp"
#include "./euler/ode_implicit_euler_jacobian_standard_policy.hpp"

#include "./base/ode_implicit_bdf2_residual_policy_base.hpp"
#include "./base/ode_implicit_bdf2_jacobian_policy_base.hpp"
#include "./bdf/ode_implicit_bdf2_residual_standard_policy.hpp"
#include "./bdf/ode_implicit_bdf2_jacobian_standard_policy.hpp"

#include "./base/ode_implicit_adams_moulton1_residual_policy_base.hpp"
#include "./base/ode_implicit_adams_moulton1_jacobian_policy_base.hpp"
#include "./adams_moulton/ode_implicit_adams_moulton1_residual_standard_policy.hpp"
#include "./adams_moulton/ode_implicit_adams_moulton1_jacobian_standard_policy.hpp"

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
	  typename ...Args>
struct isLegitimateImplicitEulerResidualPolicy<
  policy_t<state_type,residual_type,model_type,time_type,Args...>,
  typename std::enable_if<
    core::meta::publiclyInheritsFrom<
      policy_t<state_type,residual_type,model_type,time_type,Args...>,
      ode::policy::implicitEulerResidualPolicyBase<policy_t,state_type,
						   residual_type,model_type,
						   time_type, Args...
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
	  typename ...Args>
struct isLegitimateImplicitEulerJacobianPolicy<
  policy_t<state_type,residual_type,model_type,time_type,Args...>,
  typename std::enable_if<
    core::meta::publiclyInheritsFrom<
      policy_t<state_type,residual_type,model_type,time_type,Args...>,
      ode::policy::implicitEulerJacobianPolicyBase<policy_t,state_type,
						   residual_type,model_type,
						   time_type, Args...
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
	  typename time_type>
struct isImplicitEulerResidualStandardPolicy<
  policy_t<state_type,residual_type,model_type,time_type>,
  typename std::enable_if<
    std::is_same<policy_t<state_type,residual_type,model_type,time_type>,
		 ode::policy::implicitEulerStandardResidual<
		   state_type, residual_type,model_type, time_type>
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
	  typename time_type>
struct isImplicitEulerJacobianStandardPolicy<
  policy_t<state_type,jacobian_type,model_type,time_type>,
  typename std::enable_if<
    std::is_same<policy_t<state_type,jacobian_type,model_type,time_type>,
		 ode::policy::implicitEulerStandardJacobian<
		   state_type, jacobian_type,model_type, time_type>
		 >::value
			  >::type
  > : std::true_type{};
//----------------------------------------------------------------




//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//-----------------------------------------------------------
// IMPLICIT BDF2 
//-----------------------------------------------------------
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

template<typename policy_t, typename enable = void>
struct isLegitimateImplicitBDF2ResidualPolicy : std::false_type{};
  
template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type,
	  typename ...Args>
struct isLegitimateImplicitBDF2ResidualPolicy<
  policy_t<state_type,residual_type,model_type,time_type,Args...>,
  typename std::enable_if<
    core::meta::publiclyInheritsFrom<
      policy_t<state_type,residual_type,model_type,time_type,Args...>,
      ode::policy::implicitBDF2ResidualPolicyBase<policy_t,state_type,
						   residual_type,model_type,
						   time_type, Args...
						   >
                                    >::value
                         >::type
  > : std::true_type{};
//-----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isLegitimateImplicitBDF2JacobianPolicy : std::false_type{};
  
template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type,
	  typename ...Args>
struct isLegitimateImplicitBDF2JacobianPolicy<
  policy_t<state_type,residual_type,model_type,time_type,Args...>,
  typename std::enable_if<
    core::meta::publiclyInheritsFrom<
      policy_t<state_type,residual_type,model_type,time_type,Args...>,
      ode::policy::implicitBDF2JacobianPolicyBase<policy_t,state_type,
						   residual_type,model_type,
						   time_type, Args...
						   >
                                    >::value
                         >::type
  > : std::true_type{};
//-----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isImplicitBDF2ResidualStandardPolicy
  : std::false_type{};

template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type>
struct isImplicitBDF2ResidualStandardPolicy<
  policy_t<state_type,residual_type,model_type,time_type>,
  typename std::enable_if<
    std::is_same<policy_t<state_type,residual_type,model_type,time_type>,
		 ode::policy::implicitBDF2StandardResidual<
		   state_type, residual_type,model_type, time_type>
		 >::value
			  >::type
  > : std::true_type{};
//----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isImplicitBDF2JacobianStandardPolicy
  : std::false_type{};

template <template <typename...> class policy_t,
	  typename state_type,
	  typename jacobian_type,
	  typename model_type,
	  typename time_type>
struct isImplicitBDF2JacobianStandardPolicy<
  policy_t<state_type,jacobian_type,model_type,time_type>,
  typename std::enable_if<
    std::is_same<policy_t<state_type,jacobian_type,model_type,time_type>,
		 ode::policy::implicitBDF2StandardJacobian<
		   state_type, jacobian_type,model_type, time_type>
		 >::value
			  >::type
  > : std::true_type{};
//----------------------------------------------------------------




//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//-----------------------------------------------------------
// IMPLICIT ADAMS-MOULTON1
//-----------------------------------------------------------
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

template<typename policy_t, typename enable = void>
struct isLegitimateImplicitAdamsMoulton1ResidualPolicy : std::false_type{};
  
template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type,
	  typename ...Args>
struct isLegitimateImplicitAdamsMoulton1ResidualPolicy<
  policy_t<state_type,residual_type,model_type,time_type,Args...>,
  typename std::enable_if<
    core::meta::publiclyInheritsFrom<
      policy_t<state_type,residual_type,model_type,time_type,Args...>,
      ode::policy::implicitAdamsMoulton1ResidualPolicyBase<policy_t,state_type,
						   residual_type,model_type,
						   time_type, Args...
						   >
                                    >::value
                         >::type
  > : std::true_type{};
//-----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isLegitimateImplicitAdamsMoulton1JacobianPolicy : std::false_type{};
  
template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type,
	  typename ...Args>
struct isLegitimateImplicitAdamsMoulton1JacobianPolicy<
  policy_t<state_type,residual_type,model_type,time_type,Args...>,
  typename std::enable_if<
    core::meta::publiclyInheritsFrom<
      policy_t<state_type,residual_type,model_type,time_type,Args...>,
      ode::policy::implicitAdamsMoulton1JacobianPolicyBase<policy_t,state_type,
						   residual_type,model_type,
						   time_type, Args...
						   >
                                    >::value
                         >::type
  > : std::true_type{};
//-----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isImplicitAdamsMoulton1ResidualStandardPolicy
  : std::false_type{};

template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type>
struct isImplicitAdamsMoulton1ResidualStandardPolicy<
  policy_t<state_type,residual_type,model_type,time_type>,
  typename std::enable_if<
    std::is_same<policy_t<state_type,residual_type,model_type,time_type>,
		 ode::policy::implicitAdamsMoulton1StandardResidual<
		   state_type, residual_type,model_type, time_type>
		 >::value
			  >::type
  > : std::true_type{};
//----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isImplicitAdamsMoulton1JacobianStandardPolicy
  : std::false_type{};

template <template <typename...> class policy_t,
	  typename state_type,
	  typename jacobian_type,
	  typename model_type,
	  typename time_type>
struct isImplicitAdamsMoulton1JacobianStandardPolicy<
  policy_t<state_type,jacobian_type,model_type,time_type>,
  typename std::enable_if<
    std::is_same<policy_t<state_type,jacobian_type,model_type,time_type>,
		 ode::policy::implicitAdamsMoulton1StandardJacobian<
		   state_type, jacobian_type,model_type, time_type>
		 >::value
			  >::type
  > : std::true_type{};
//----------------------------------------------------------------
  

  
} // namespace meta
} // namespace core
#endif
