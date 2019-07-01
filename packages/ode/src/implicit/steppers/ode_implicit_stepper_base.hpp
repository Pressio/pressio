
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_BASE_IMPLICIT_STEPPER_BASE_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_BASE_IMPLICIT_STEPPER_BASE_HPP_

#include "ode_implicit_stepper_traits.hpp"
#include "../policies/meta/ode_is_implicit_jacobian_standard_policy.hpp"
#include "../policies/meta/ode_is_implicit_residual_standard_policy.hpp"
#include "../policies/meta/ode_is_legitimate_implicit_jacobian_policy.hpp"
#include "../policies/meta/ode_is_legitimate_implicit_residual_policy.hpp"
#include "../../ode_storage.hpp"
#include "../../ode_system_wrapper.hpp"
//#include "../ode_implicit_aux_data.hpp"

namespace rompp{ namespace ode{

/*
 * (1) constructors here should be private but we need
 * them public to enable interfacing with pybind11
 */

template<typename concrete_stepper_type, int nAuxStates>
class ImplicitStepperBase
//: private utils::details::CrtpBase<ImplicitStepperBase<concrete_stepper_type, nAuxStates>>
{
  using traits		  = typename details::traits<concrete_stepper_type>;
  using sc_t		  = typename traits::scalar_t;
  using state_t		  = typename traits::state_t;
  using residual_t	  = typename traits::residual_t;
  using jacobian_t	  = typename traits::jacobian_t;
  using standard_res_policy_t = typename traits::standard_res_policy_t;
  using standard_jac_policy_t = typename traits::standard_jac_policy_t;
  using residual_pol_t = typename traits::residual_policy_t;
  using jacobian_pol_t = typename traits::jacobian_policy_t;
  using model_t		  = typename traits::model_t;
  using system_wrapper_t = impl::OdeSystemWrapper<model_t>;

  //do checking here that things are as supposed
  static_assert( meta::is_legitimate_implicit_state_type<state_t>::value,
       "OOPS: STATE_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");
  static_assert( meta::is_legitimate_implicit_residual_type<residual_t>::value,
       "OOPS: RESIDUAL_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");
  static_assert( meta::is_legitimate_jacobian_type<jacobian_t>::value,
       "OOPS: JACOBIAN_TYPE IN SELECTED IMPLICIT STEPPER IS NOT VALID");

protected:
  // procted because these are accessed only by children classes
  sc_t t_  = {};
  sc_t dt_ = {};
  system_wrapper_t sys_;
  impl::OdeStorage<state_t, nAuxStates> stateAuxStorage_;

  // conditionally set the type of the object knowing how to compute residual
  // if we have a standard policy, then it takes a copy
  // if we have a user-defined policy, we take a const & to it
  typename std::conditional<
    mpl::is_same<standard_res_policy_t, residual_pol_t>::value,
    const residual_pol_t,
    const residual_pol_t &
    >::type residual_obj_;

  // conditionally set the type of the object knowing how to compute jacobian
  // if we have a standard policy, then it takes a copy
  // if we have a user-defined policy, we take a const & to it
  typename std::conditional<
    mpl::is_same<standard_jac_policy_t, jacobian_pol_t>::value,
    const jacobian_pol_t,
    const jacobian_pol_t &
    >::type jacobian_obj_;

public:
  decltype(traits::order_value) order() const{
    return traits::order_value;
  }

  void residual(const state_t & y,
  		residual_t & R) const{
    this->residual_obj_.template operator()<
      traits::enum_id,
      traits::steps
      >(y, R, stateAuxStorage_.data_, sys_.get(), this->t_, this->dt_);
  }

  residual_t residual(const state_t & y) const{
    std::cout << " residual_impl_st_base" << std::endl;
    return this->residual_obj_.template operator()<
      traits::enum_id,
      traits::steps
      >(y, stateAuxStorage_.data_, sys_.get(), this->t_, this->dt_);
  }

  void jacobian(const state_t & y,
  		jacobian_t & J) const{
    this->jacobian_obj_.template operator()<
      traits::enum_id
      >(y, J, sys_.get(), this->t_, this->dt_);
  }

  jacobian_t jacobian(const state_t & y) const{
    return this->jacobian_obj_.template operator()<
      traits::enum_id
      >(y, sys_.get(), this->t_, this->dt_);
  }

public:
  ImplicitStepperBase() = delete;
  ~ImplicitStepperBase() = default;

  ImplicitStepperBase(const state_t & y0,
		      const model_t & model,
		      const residual_pol_t & resPolicyObj,
		      const jacobian_pol_t & jacPolicyObj)
    : stateAuxStorage_{y0},
      sys_{model},
      residual_obj_{resPolicyObj},
      jacobian_obj_{jacPolicyObj}{}

  // cstr for standard residual and jacob policies
  template <
    typename T1 = standard_res_policy_t,
    typename T2 = standard_jac_policy_t,
    ::rompp::mpl::enable_if_t<
      mpl::is_same<T1, residual_pol_t>::value and
      mpl::is_same<T2, jacobian_pol_t>::value
      > * = nullptr
    >
  ImplicitStepperBase(const state_t & y0,
  		      const model_t & model)
    : stateAuxStorage_{y0},
      sys_{model},
      residual_obj_{},
      jacobian_obj_{}{
	std::cout << "base stepper cnstr" << std::endl;
      }

  // cstr for standard jacob policies
  template <
    typename T2 = standard_jac_policy_t,
    ::rompp::mpl::enable_if_t<
      mpl::is_same<T2, jacobian_pol_t>::value
      > * = nullptr
    >
  ImplicitStepperBase(const state_t & y0,
  		      const model_t & model,
  		      const residual_pol_t & resPolicyObj)
    : stateAuxStorage_{y0},
      sys_{model},
      residual_obj_{resPolicyObj},
      jacobian_obj_{}{}

  // /* workaround for nvcc issue with templates, see https://devtalk.nvidia.com/default/topic/1037721/nvcc-compilation-error-with-template-parameter-as-a-friend-within-a-namespace/ */
  // template<typename DummyType> struct dummy{using type = DummyType;};
  // friend typename dummy<concrete_stepper_type>::type;
  // friend utils::details::CrtpBase<ImplicitStepperBase<concrete_stepper_type, nAuxStates>>;

};//end class

}}//end namespace rompp::ode
#endif
