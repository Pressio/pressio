
#ifndef ODE_STORAGE_HPP_
#define ODE_STORAGE_HPP_

#include "ode_ConfigDefs.hpp"

namespace rompp{ namespace ode{ namespace impl{

template<typename state_type, typename rhs_type,
	 int numAuxStates, int numAuxRHS = 0>
struct OdeStorage;
//--------------------------------------------------

template<typename state_type, typename rhs_type>
struct OdeStorage<state_type, rhs_type, 1, 0>{
  OdeStorage(state_type const & y)
    : auxStates_{{y}}{}

  ~OdeStorage() = default;

  std::array<state_type, 1> auxStates_;
};
//--------------------------------------------------

template<typename state_type, typename rhs_type>
struct OdeStorage<state_type, rhs_type, 2, 0>{

  OdeStorage(state_type const & y)
    : auxStates_{{y, y}}{}

  ~OdeStorage() = default;

  std::array<state_type, 2> auxStates_;
};
//--------------------------------------------------

template<typename state_type, typename rhs_type>
struct OdeStorage<state_type, rhs_type, 1, 4>{
  using rhs_wrapped_t = typename core::details::traits<rhs_type>::wrapped_t;

  OdeStorage(state_type const & y,
	     rhs_type const & r)
    : auxStates_{{y}},
      auxRHS_{{r, r, r, r}}
  {}

  ~OdeStorage() = default;

  std::array<state_type, 1> auxStates_;
  std::array<rhs_type, 4> auxRHS_;
};
//--------------------------------------------------

template<typename state_type, typename rhs_type>
struct OdeStorage<state_type, rhs_type, 0, 1>{
  using rhs_wrapped_t = typename core::details::traits<rhs_type>::wrapped_t;

  OdeStorage(rhs_type const & r)
    : auxRHS_{{r}}{}

  ~OdeStorage() = default;

  std::array<rhs_type, 1> auxRHS_;
};
//--------------------------------------------------

}}}//end namespace rompp::ode::impl
#endif
