
#ifndef ODE_STORAGE_HPP_
#define ODE_STORAGE_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode{
namespace impl{

template<typename state_type,
	 typename residual_type,
	 int numAuxStates,
	 int numAuxRHS = 0>
class odeStorage;

//--------------------------------------------------

template<typename state_type, typename residual_type>
class odeStorage<state_type, residual_type, 1>{
public:
  template<typename... Args>
  odeStorage(Args&&... rest)
    : auxStates_(std::forward<Args>(rest)...){}

protected:
  std::array<state_type, 1> auxStates_;
};

//--------------------------------------------------

template<typename state_type, typename residual_type>
class odeStorage<state_type, residual_type, 1, 4>{
public:
  template<typename... Args>
  odeStorage(Args&&... rest)
    : auxStates_(std::forward<Args>(rest)...),
      auxRHS_(std::forward<Args>(rest)...,
	      std::forward<Args>(rest)...,
	      std::forward<Args>(rest)...,
	      std::forward<Args>(rest)...){}

protected:
  std::array<state_type, 1> auxStates_;
  std::array<residual_type, 4> auxRHS_;
};

//--------------------------------------------------
  
template<typename state_type, typename residual_type>
class odeStorage<state_type, residual_type, 0, 1>{
public:
  template<typename... Args>
  odeStorage(Args&&... rest)
    : auxRHS_(std::forward<Args>(rest)...){}

protected:
  std::array<residual_type, 1> auxRHS_;
};

//--------------------------------------------------


}//end namespace impl
}//end namespace ode  
  
#endif
  
