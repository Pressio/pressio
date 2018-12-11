
#ifndef ODE_POLICIES_BASE_JACOBIAN_POLICY_BASE_HPP_
#define ODE_POLICIES_BASE_JACOBIAN_POLICY_BASE_HPP_

#include "../../../ode_ConfigDefs.hpp"

namespace rompp{ namespace ode{ namespace policy{

template <typename derived_t>
class JacobianPolicyBase
  : private core::details::CrtpBase<
  JacobianPolicyBase<derived_t>>{
public:

  template <typename state_type,
	    typename jacobian_type,
	    typename model_type,
	    typename scalar_type>
  void operator()(const state_type & y,
		  jacobian_type & J,
		  const model_type & model,
		  scalar_type t,
		  scalar_type dt) const
  {
    this->underlying()(y, J, model, t, dt);
  }
  //----------------------------------------

  template <typename state_type,
	    typename model_type,
	    typename scalar_type>
  auto operator()(const state_type & y,
		  const model_type & model,
		  scalar_type t,
		  scalar_type dt) const
   -> decltype(this->underlying()(y, model, t, dt))
  {
    return this->underlying()(y, model, t, dt);
  }
  //------------------------------------------

private:
  friend derived_t;
  friend core::details::CrtpBase<
    JacobianPolicyBase<derived_t>>;

  JacobianPolicyBase() = default;
  ~JacobianPolicyBase() = default;

};//end class

}//end namespace polices
}//end namespace ode
}//end namespace rompp
#endif
