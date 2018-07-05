
#ifndef ODE_ADVANCE_FULL_STATE_POLICY_BASE_HPP_
#define ODE_ADVANCE_FULL_STATE_POLICY_BASE_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace policy{


template <template <typename...> class derived_type,
	  typename ... Args>
class advanceFullStatePolicyBase
{
public:
  static constexpr bool advanceIncrement = false;
  static constexpr bool advanceFull = !advanceIncrement;

private:
  using derived_t = derived_type<Args...>;
  friend derived_t; 

  advanceFullStatePolicyBase() = default;
  ~advanceFullStatePolicyBase() = default;

};//end class

}//end namespace polices
}//end namespace ode  
#endif 
