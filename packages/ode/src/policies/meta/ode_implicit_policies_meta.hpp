
#ifndef ODE_IMPLICIT_POLICIES_META_HPP_
#define ODE_IMPLICIT_POLICIES_META_HPP_

#include "../base/ode_residual_policy_base.hpp"
#include "../base/ode_jacobian_policy_base.hpp"

#include "../standard/ode_residual_standard_policy.hpp"
#include "../standard/ode_jacobian_standard_policy.hpp"

namespace ode{
namespace meta {


//////////////////////////////////////////////////////////////////
//
// IMPLICIT EULER 
//
//////////////////////////////////////////////////////////////////

template<typename policy_t, typename enable = void>
struct isLegitimateImplicitEulerResidualPolicy 
: isLegitimateResidualPolicy<policy_t>{};

//-----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isLegitimateImplicitEulerJacobianPolicy 
: isLegitimateJacobianPolicy<policy_t>{};

//-----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isImplicitEulerResidualStandardPolicy 
: isResidualStandardPolicy<policy_t>{};

//----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isImplicitEulerJacobianStandardPolicy 
: isJacobianStandardPolicy<policy_t>{};

//----------------------------------------------------------------


//////////////////////////////////////////////////////////////////
//
// IMPLICIT BDF2
//
//////////////////////////////////////////////////////////////////

template<typename policy_t, typename enable = void>
struct isLegitimateImplicitBDF2ResidualPolicy 
: isLegitimateResidualPolicy<policy_t>{};

//-----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isLegitimateImplicitBDF2JacobianPolicy 
: isLegitimateJacobianPolicy<policy_t>{};

//-----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isImplicitBDF2ResidualStandardPolicy 
: isResidualStandardPolicy<policy_t>{};

//----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isImplicitBDF2JacobianStandardPolicy 
: isJacobianStandardPolicy<policy_t>{};

//----------------------------------------------------------------


//////////////////////////////////////////////////////////////////
//
// IMPLICIT BDF3
//
//////////////////////////////////////////////////////////////////

template<typename policy_t, typename enable = void>
struct isLegitimateImplicitBDF3ResidualPolicy 
: isLegitimateResidualPolicy<policy_t>{};

//-----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isLegitimateImplicitBDF3JacobianPolicy 
: isLegitimateJacobianPolicy<policy_t>{};

//-----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isImplicitBDF3ResidualStandardPolicy 
: isResidualStandardPolicy<policy_t>{};

//----------------------------------------------------------------

template<typename policy_t, typename enable = void>
struct isImplicitBDF3JacobianStandardPolicy 
: isJacobianStandardPolicy<policy_t>{};

//----------------------------------------------------------------

  
} // namespace meta
} // namespace core
#endif
