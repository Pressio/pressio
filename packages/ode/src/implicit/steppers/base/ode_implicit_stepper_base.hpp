
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_BASE_IMPLICIT_STEPPER_BASE_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_BASE_IMPLICIT_STEPPER_BASE_HPP_

#include "../ode_implicit_stepper_traits.hpp"
#include "../../../ode_basic_meta.hpp"
#include "../../policies/meta/ode_is_implicit_jacobian_standard_policy.hpp"
#include "../../policies/meta/ode_is_implicit_residual_standard_policy.hpp"
#include "../../policies/meta/ode_is_legitimate_implicit_jacobian_policy.hpp"
#include "../../policies/meta/ode_is_legitimate_implicit_residual_policy.hpp"
#include "../../../ode_storage.hpp"
#include "../../../ode_aux_data.hpp"

namespace rompp{ namespace ode{

template<typename concrete_stepper_type>
class ImplicitStepperBase
  : private core::details::CrtpBase<ImplicitStepperBase<concrete_stepper_type>>
{
  using traits		  = typename details::traits<concrete_stepper_type>;
  using sc_t		  = typename traits::scalar_t;
  using state_t		  = typename traits::state_t;
  using residual_t	  = typename traits::residual_t;
  using jacobian_t	  = typename traits::jacobian_t;
  using residual_policy_t = typename traits::residual_policy_t;
  using jacobian_policy_t = typename traits::jacobian_policy_t;

  //do checking here that things are as supposed
  static_assert( meta::is_legitimate_implicit_state_type<state_t>::value,
       "OOPS: STATE_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");
  static_assert( meta::is_legitimate_implicit_residual_type<residual_t>::value,
       "OOPS: RESIDUAL_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");
  static_assert( meta::is_legitimate_jacobian_type<jacobian_t>::value,
       "OOPS: JACOBIAN_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");

public:
  typename traits::order_t order() const{
    return traits::order_value;
  }

  void residual(const state_t & y, residual_t & R)const{
    this->underlying().residualImpl(y, R);
  }

  void jacobian(const state_t & y, jacobian_t & J)const{
    this->underlying().jacobianImpl(y, J);
  }

  residual_t residual(const state_t & y)const{
    return this->underlying().residualImpl(y);
  }

  jacobian_t jacobian(const state_t & y)const{
    return this->underlying().jacobianImpl(y);
  }

private:
  ImplicitStepperBase() = default;
  ~ImplicitStepperBase() = default;

  friend concrete_stepper_type;
  friend core::details::CrtpBase<ImplicitStepperBase<concrete_stepper_type>>;

};//end class

}}//end namespace rompp::ode
#endif
