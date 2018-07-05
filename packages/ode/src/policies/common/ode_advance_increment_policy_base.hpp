
#ifndef ODE_ADVANCE_INCREMENT_POLICY_BASE_HPP_
#define ODE_ADVANCE_INCREMENT_POLICY_BASE_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace policy{

template <template <typename...> class derived_type,
	  typename state_type,
	  typename ... Args>
class advanceIncrementPolicyBase
{
public:
  static constexpr bool advanceIncrement = true;
  static constexpr bool advanceFull = !advanceIncrement;

protected:
  state_type const * y0ptr_;
  state_type yFull_;  

private:
  using derived_t = derived_type<state_type, Args...>;
  friend derived_t; 

  advanceIncrementPolicyBase(const state_type & y0)
    : y0ptr_(&y0), yFull_(y0){}

  ~advanceIncrementPolicyBase() = default;  

};//end class

}//end namespace polices
}//end namespace ode  
#endif 
