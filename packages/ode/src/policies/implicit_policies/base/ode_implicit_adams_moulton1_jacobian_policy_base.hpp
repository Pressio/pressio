
#ifndef ODE_IMPLICIT_ADAMS_MOULTON1_JACOBIAN_POLICY_BASE_HPP_
#define ODE_IMPLICIT_ADAMS_MOULTON1_JACOBIAN_POLICY_BASE_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace policy{
  
template <template <typename...> class derived_type,
	  typename state_type,
	  typename jacobian_type,
	  typename model_type,
	  typename time_type,
	  typename ... Args>
class implicitAdamsMoulton1JacobianPolicyBase
{
public:
  void compute(const state_type & y,
	       jacobian_type & J,
	       model_type & model,
	       time_type t,
	       time_type dt){
    this->underlying().computeImpl(y, J, model, t, dt);
  } 

private:
  using derived_t = derived_type<state_type,jacobian_type,
				 model_type, time_type, Args...>;
  friend derived_t; 
  implicitAdamsMoulton1JacobianPolicyBase() = default;
  ~implicitAdamsMoulton1JacobianPolicyBase() = default;
  
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
