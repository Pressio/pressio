
#ifndef ODE_EXPLICIT_ANY_RUNGEKUTTA4_STEPPER_IMPL_HPP_
#define ODE_EXPLICIT_ANY_RUNGEKUTTA4_STEPPER_IMPL_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace impl{


template<typename state_type,
	 typename residual_type,
	 int stages>
class rungeKuttaStorage;

  
template<typename state_type, typename residual_type>
class rungeKuttaStorage<state_type, residual_type, 4>
{
public:
  template<typename... Args>
  rungeKuttaStorage(Args&&... rest)
    : rhs_(std::forward<Args>(rest)...,
	   std::forward<Args>(rest)...,
	   std::forward<Args>(rest)...,
	   std::forward<Args>(rest)...)
  {}

protected:
  std::array<residual_type, 4> rhs_;
};


}//end namespace impl
}//end namespace ode  
  
#endif
  
