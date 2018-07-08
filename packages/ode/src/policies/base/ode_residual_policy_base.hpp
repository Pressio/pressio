
#ifndef ODE_RESIDUAL_POLICY_BASE_HPP_
#define ODE_RESIDUAL_POLICY_BASE_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace policy{
    
template <template <typename...> class derived_type,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type,
	  typename sizer_type, 
	  typename ... Args>
class residualPolicyBase
{
public:
  void compute(const state_type & y,
	       residual_type & R,
	       model_type & model,
	       time_type t){
    this->underlying().computeImpl(y, R, model, t);
  }

  void weightTimeDiscreteResidual(const state_type & y,
         residual_type & R,
         model_type & model, 
         time_type t){
    this->underlying().weightTimeDiscreteResidualImpl(y, R, model, t);
  }

private:
  using derived_t = derived_type<state_type,
				 residual_type,
				 model_type,
				 time_type,
				 sizer_type,
				 Args...>;
  friend derived_t;

  residualPolicyBase() = default;
  ~residualPolicyBase() = default;
  
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
