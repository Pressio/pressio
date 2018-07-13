
#ifndef ODE_JACOBIAN_POLICY_BASE_HPP_
#define ODE_JACOBIAN_POLICY_BASE_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace policy{
  
template <typename derived_t>
class jacobianPolicyBase
{
public:

  template <typename state_type, typename jacobian_type,
	    typename model_type, typename time_type>
  void compute(const state_type & y, 
	       jacobian_type & J,
	       model_type & model, 
	       time_type t,
	       time_type dt)
  {
    this->underlying().computeImpl(y, J, model, t, dt);
  } 
  
private:
  friend derived_t;

  jacobianPolicyBase() = default;
  ~jacobianPolicyBase() = default;
  
  derived_t & underlying(){
    return static_cast<derived_t& >( *this );
  }
  derived_t const & underlying() const{
    return static_cast<derived_t const & >( *this );
  }
};


}//end namespace polices
}//end namespace ode  
#endif 
