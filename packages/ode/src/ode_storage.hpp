
#ifndef ODE_STORAGE_HPP_
#define ODE_STORAGE_HPP_

#include "ode_ConfigDefs.hpp"

namespace rompp{ namespace ode{ namespace impl{

template<typename state_type, typename rhs_type,
	 int numAuxStates, int numAuxRHS = 0>
class OdeStorage;
//--------------------------------------------------

template<typename state_type, typename rhs_type>
class OdeStorage<state_type, rhs_type, 1, 0>{
public:
  OdeStorage(state_type const & y)
    : auxStates_{{y}}{}

  ~OdeStorage() = default;

protected:
  std::array<state_type, 1> auxStates_;
};
//--------------------------------------------------

template<typename state_type, typename rhs_type>
class OdeStorage<state_type, rhs_type, 2, 0>{
public:
  OdeStorage(state_type const & y)
    : auxStates_{{y, y}}{}

  ~OdeStorage() = default;

protected:
  std::array<state_type, 2> auxStates_;
};
//--------------------------------------------------


template<typename state_type, typename rhs_type>
class OdeStorage<state_type, rhs_type, 1, 4>{
  using rhs_wrapped_t = typename core::details::traits<rhs_type>::wrapped_t;

public:
  OdeStorage(state_type const & y,
	     rhs_type const & r)
    : auxStates_{{y}},
      auxRHS_{{r, r, r, r}}
  {}

  ~OdeStorage() = default;

protected:
  std::array<state_type, 1> auxStates_;
  std::array<rhs_type, 4> auxRHS_;
};
//--------------------------------------------------


template<typename state_type, typename rhs_type>
class OdeStorage<state_type, rhs_type, 0, 1>{
  using rhs_wrapped_t = typename core::details::traits<rhs_type>::wrapped_t;

public:
  OdeStorage(rhs_type const & r) : auxRHS_{r}{}

  // OdeStorage(rhs_wrapped_t const & rNative)
  //   : auxRHS_{rhs_type{rNative}}{}

  ~OdeStorage() = default;

protected:
  std::array<rhs_type, 1> auxRHS_;
};
//--------------------------------------------------

}}}//end namespace rompp::ode::impl
#endif
