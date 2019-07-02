
#ifndef ODE_POLICIES_STANDARD_EXPLICIT_VELOCITY_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_EXPLICIT_VELOCITY_STANDARD_POLICY_HPP_

#include "../../ode_fwd.hpp"
#include "ode_explicit_velocity_policy_base.hpp"

#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#endif

namespace pressio{ namespace ode{ namespace policy{

/*
 * state_type = velocity_type
 * both are wrappers from containers
 */
template<
  typename state_type,
  typename system_type
  >
class ExplicitVelocityStandardPolicy<
  state_type, system_type, state_type,
  mpl::enable_if_t<
    containers::meta::is_vector_wrapper<state_type>::value
#ifdef HAVE_PYBIND11
    and mpl::not_same<system_type, pybind11::object >::value
#endif
    >
  >
  : public ExplicitVelocityPolicyBase<
  ExplicitVelocityStandardPolicy<state_type, system_type> >{

  using base_t = ExplicitVelocityPolicyBase<
    ExplicitVelocityStandardPolicy<state_type, system_type>>;
  friend base_t;

public:
  ExplicitVelocityStandardPolicy() = default;
  ~ExplicitVelocityStandardPolicy() = default;

  template < typename scalar_type >
  void operator()(const state_type & y,
		  state_type & R,
		  const system_type & model,
		  scalar_type t) const{
    model.velocity(*y.data(), *R.data(), t);
  }

  template < typename scalar_type >
  state_type operator()(const state_type & y,
			const system_type & model,
			scalar_type t) const{
    return state_type(model.velocity(*y.data(), t));
  }
};//end class



/*
 * state_type = velocity_type
 * both are pybind11::array_t
 */
#ifdef HAVE_PYBIND11
template<
  typename state_type,
  typename system_type
  >
class ExplicitVelocityStandardPolicy<
  state_type, system_type, state_type,
  mpl::enable_if_t<
    mpl::is_same<system_type, pybind11::object >::value and
    containers::meta::is_cstyle_array_pybind11<state_type>::value
    >
  >
  : public ExplicitVelocityPolicyBase<
  ExplicitVelocityStandardPolicy<state_type, system_type> >{

  using base_t = ExplicitVelocityPolicyBase<
    ExplicitVelocityStandardPolicy<state_type, system_type>>;
  friend base_t;

public:
  ExplicitVelocityStandardPolicy() = default;
  ~ExplicitVelocityStandardPolicy() = default;

  template <typename scalar_type>
  void operator()(const state_type & y,
		  state_type & R,
		  const system_type & model,
		  scalar_type t) const{
    //printf("C++ R address: %p\n", R.data());
    model.attr("velocity2")(y, R, t);
  }

  template <typename scalar_type>
  state_type operator()(const state_type & y,
  			const system_type & model,
  			scalar_type t) const{
    return model.attr("velocity1")(y, t);
  }

};//end class
#endif


}}}//end namespace pressio::ode::policy
#endif
