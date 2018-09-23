
#ifndef ODE_POLICIES_META_EXPLICIT_RUNGE_KUTTA4_POLICIES_META_HPP_
#define ODE_POLICIES_META_EXPLICIT_RUNGE_KUTTA4_POLICIES_META_HPP_

#include "ode_explicit_policies_meta.hpp"
#include "../standard/ode_explicit_residual_standard_policy.hpp"

namespace rompp{
namespace ode{
namespace meta {

//-----------------------------------------------------------------
// METAF TO CHECK RESIDUAL POLICY IS STANDARD
//-----------------------------------------------------------------
template<typename policy_t, typename enable = void>
struct is_explicit_runge_kutta4_residual_standard_policy
  : std::false_type{};

template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type>
struct is_explicit_runge_kutta4_residual_standard_policy<
  policy_t<state_type, residual_type, model_type>,
  typename std::enable_if<
    std::is_same<policy_t<state_type, residual_type, model_type>,
		 ode::policy::ExplicitResidualStandardPolicy<
		   state_type, residual_type, model_type>
		 >::value
    >::type
  > : std::true_type{};
//----------------------------------------------------------------

} // namespace meta
} // namespace core
}//end namespace rompp
#endif
