
#ifndef ODE_STORAGE_HPP_
#define ODE_STORAGE_HPP_

#include "ode_ConfigDefs.hpp"

namespace rompp{
namespace ode{
namespace impl{

template<typename state_type, typename rhs_type,
	 int numAuxStates, int numAuxRHS = 0>
class OdeStorage;

//--------------------------------------------------

template<typename state_type, typename rhs_type>
class OdeStorage<state_type, rhs_type, 1>{
public:
  OdeStorage(state_type const & y0)
    : auxStates_{y0}{}
  ~OdeStorage() = default;

protected:
  std::array<state_type, 1> auxStates_;
};
//--------------------------------------------------

template<typename state_type, typename rhs_type>
class OdeStorage<state_type, rhs_type, 1, 4>{
public:
  OdeStorage(state_type const & y0,
	     rhs_type const & r0)
    : auxStates_{y0},
      auxRHS_{r0, r0, r0, r0}{}
  ~OdeStorage() = default;

protected:
  std::array<state_type, 1> auxStates_;
  std::array<rhs_type, 4> auxRHS_;
};
//--------------------------------------------------
  
template<typename state_type, typename rhs_type>
class OdeStorage<state_type, rhs_type, 0, 1>{
public:
  OdeStorage(rhs_type const & r0)
    :  auxRHS_{r0}{}
  ~OdeStorage() = default;

protected:
  std::array<rhs_type, 1> auxRHS_;
};
//--------------------------------------------------

}//end namespace impl
}//end namespace ode  
  
}//end namespace rompp
#endif
  
