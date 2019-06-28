
#ifndef ODE_POLICIES_STANDARD_EXPLICIT_RESIDUAL_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_EXPLICIT_RESIDUAL_STANDARD_POLICY_HPP_

#include "../../ode_forward_declarations.hpp"
#include "ode_explicit_residual_policy_base.hpp"

#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#endif

namespace rompp{ namespace ode{ namespace policy{

/*
 * state_type = residual_type
 * both are wrappers from algebra
 */
template<
  typename state_type,
  typename model_type
  >
class ExplicitResidualStandardPolicy<
  state_type,model_type, state_type,
  mpl::enable_if_t<
    algebra::meta::is_algebra_vector_wrapper<state_type>::value
#ifdef HAVE_PYBIND11
    and mpl::not_same<model_type, pybind11::object >::value
#endif
    >
  >
  : public ExplicitResidualPolicyBase<
  ExplicitResidualStandardPolicy<state_type, model_type> >{

  using base_t = ExplicitResidualPolicyBase<
    ExplicitResidualStandardPolicy<state_type, model_type>>;
  friend base_t;

public:
  ExplicitResidualStandardPolicy() = default;
  ~ExplicitResidualStandardPolicy() = default;

  template < typename scalar_type >
  void operator()(const state_type & y,
		  state_type & R,
		  const model_type & model,
		  scalar_type t) const{
    model.residual(*y.data(), *R.data(), t);
  }

  template < typename scalar_type >
  state_type operator()(const state_type & y,
			const model_type & model,
			scalar_type t) const{
    return state_type(model.residual(*y.data(), t));
  }
};//end class



/*
 * state_type = residual_type
 * both are pybind11::array_t
 */
#ifdef HAVE_PYBIND11
template<
  typename state_type,
  typename model_type
  >
class ExplicitResidualStandardPolicy<
  state_type,model_type, state_type,
  mpl::enable_if_t<
    mpl::is_same<model_type, pybind11::object >::value and
    algebra::meta::is_cstyle_array_pybind11<state_type>::value
    >
  >
  : public ExplicitResidualPolicyBase<
  ExplicitResidualStandardPolicy<state_type, model_type> >{

  using base_t = ExplicitResidualPolicyBase<
    ExplicitResidualStandardPolicy<state_type, model_type>>;
  friend base_t;

public:
  ExplicitResidualStandardPolicy() = default;
  ~ExplicitResidualStandardPolicy() = default;

  template <typename scalar_type>
  void operator()(const state_type & y,
		  state_type & R,
		  const model_type & model,
		  scalar_type t) const{
    //printf("C++ R address: %p\n", R.data());
    model.attr("residual2")(y, R, t);
  }

  template <typename scalar_type>
  state_type operator()(const state_type & y,
  			const model_type & model,
  			scalar_type t) const{
    return model.attr("residual1")(y, t);
  }

};//end class
#endif


}}}//end namespace rompp::ode::policy
#endif
