
#ifndef ODE_JACOBIAN_POLICY_BASE_HPP_
#define ODE_JACOBIAN_POLICY_BASE_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace policy{
  
template <template <typename...> class derived_type,
	  typename state_type,
	  typename jacobian_type,
	  typename model_type,
	  typename time_type,
    typename sizer_type,
	  typename ... Args>
class jacobianPolicyBase
{
public:
  void compute(const state_type & y, 
	       jacobian_type & J,
	       model_type & model, 
	       time_type t){
    this->underlying().computeImpl(y, J, model, t);
  } 

  void weightJacobian(const state_type & y,
		      jacobian_type & J,
		      model_type & model, 
		      time_type t){
    this->underlying().weightJacobianImpl(y, J, model, t);
  }
  
private:
  using derived_t = derived_type<state_type,
                  jacobian_type,
				          model_type, 
                  time_type, 
                  sizer_type,
                  Args...>;

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
