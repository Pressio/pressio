
#ifndef ODE_STORAGE_HPP_
#define ODE_STORAGE_HPP_

#include "ode_ConfigDefs.hpp"
#include <array>
#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#endif

namespace rompp{ namespace ode{ namespace impl{

/*
 * note that these are auxiliary objects for storing data.
 * it is fundamental that these DO not point to the same memory
 * locations of the objects passed in.
 * In other words, the y,r arguments to constructors are
 * only needed so that we copy-construct them since
 * we do not know if they have or not a default constructor, etc.
 * each type can be different so we need a general way to create
 * an object from one of the same kind.
 *
 * However, if copy-constructing implements shallow copy, then
 * we need to do something else because the storage objects
 * need to be new allocations.
 */

template<
  typename state_type,
  typename rhs_type,
  int numAuxStates,
  int numAuxRHS = 0
  >
struct OdeStorage;


//--------------------------------------------------
// num_aux_states = 1, num_aux_rhs = 0
//--------------------------------------------------
template<typename state_type, typename rhs_type>
struct OdeStorage<state_type, rhs_type, 1, 0>{

  template <
    typename _state_type = state_type
#ifdef HAVE_PYBIND11
    , mpl::enable_if_t<
      !containers::meta::is_array_pybind11<_state_type>::value
      > * = nullptr
#endif
    >
  OdeStorage(_state_type const & y)
    : auxStates_{{y}}{}

#ifdef HAVE_PYBIND11
  template <
    typename _state_type = state_type,
    mpl::enable_if_t<
      containers::meta::is_array_pybind11<_state_type>::value
      > * = nullptr
    >
  OdeStorage(_state_type const & y)
    : auxStates_{{_state_type(const_cast<_state_type &>(y).request())}}
  {}
#endif

  ~OdeStorage() = default;

  std::array<state_type, 1> auxStates_;
};



//--------------------------------------------------
// num_aux_states = 2, num_aux_rhs = 0
//--------------------------------------------------
template<typename state_type, typename rhs_type>
struct OdeStorage<state_type, rhs_type, 2, 0>{

  template <
    typename _state_type = state_type
#ifdef HAVE_PYBIND11
    , mpl::enable_if_t<
      !containers::meta::is_array_pybind11<_state_type>::value
      > * = nullptr
#endif
    >
  OdeStorage(_state_type const & y)
    : auxStates_{{y,y}}{}

#ifdef HAVE_PYBIND11
  template <
    typename _state_type = state_type,
    mpl::enable_if_t<
      containers::meta::is_array_pybind11<_state_type>::value
      > * = nullptr
    >
  OdeStorage(_state_type const & y)
    : auxStates_{{_state_type(const_cast<_state_type &>(y).request()),
		  _state_type(const_cast<_state_type &>(y).request())}}
  {}
#endif

  ~OdeStorage() = default;

  std::array<state_type, 2> auxStates_;
};



//--------------------------------------------------
// num_aux_states = 1, num_aux_rhs = 4
//--------------------------------------------------
template<typename state_type, typename rhs_type>
struct OdeStorage<state_type, rhs_type, 1, 4>{

  template <
    typename _state_type = state_type,
    typename _rhs_type = rhs_type
#ifdef HAVE_PYBIND11
    ,mpl::enable_if_t<
       !containers::meta::is_array_pybind11<_state_type>::value and
       !containers::meta::is_array_pybind11<_rhs_type>::value
       > * = nullptr
#endif
    >
  OdeStorage(_state_type const & y,
	     _rhs_type const & r)
    : auxStates_{{y}},
      auxRHS_{{r, r, r, r}}{}

#ifdef HAVE_PYBIND11
  template <
    typename _state_type = state_type,
    typename _rhs_type = rhs_type,
    mpl::enable_if_t<
      containers::meta::is_array_pybind11<_state_type>::value and
      containers::meta::is_array_pybind11<_rhs_type>::value
      > * = nullptr
    >
    OdeStorage(_state_type const & y, _rhs_type r)
      // this syntax bascially creates new object and copies values
      // because we DO not want to do shallow copies
      : auxStates_{_state_type(const_cast<_state_type &>(y).request())},
	auxRHS_{{_rhs_type(const_cast<_rhs_type &>(r).request()),
		 _rhs_type(const_cast<_rhs_type &>(r).request()),
		 _rhs_type(const_cast<_rhs_type &>(r).request()),
		 _rhs_type(const_cast<_rhs_type &>(r).request())}}
  {
    printf("OdeStorage14 y          addr: %p\n", y.data());
    printf("OdeStorage14 r          addr: %p\n", r.data());
    printf("OdeStorage14 auxRHS_[0] addr: %p\n", auxRHS_[0].data());
    printf("OdeStorage14 auxRHS_[1] addr: %p\n", auxRHS_[1].data());
    printf("OdeStorage14 auxRHS_[2] addr: %p\n", auxRHS_[2].data());
    printf("OdeStorage14 auxRHS_[3] addr: %p\n", auxRHS_[3].data());
  }
#endif

  ~OdeStorage() = default;

  std::array<state_type, 1> auxStates_;
  std::array<rhs_type, 4> auxRHS_;
};


//--------------------------------------------------
// num_aux_states = 0, num_aux_rhs = 1
//--------------------------------------------------
template<typename state_type, typename rhs_type>
struct OdeStorage<state_type, rhs_type, 0, 1>{

  template <
    typename _rhs_type = rhs_type
#ifdef HAVE_PYBIND11
    , mpl::enable_if_t<
	!containers::meta::is_array_pybind11<_rhs_type>::value
      > * = nullptr
#endif
    >
  OdeStorage(_rhs_type const & r)
    : auxRHS_{{r}}{}

#ifdef HAVE_PYBIND11
  template <
    typename _rhs_type = rhs_type,
    mpl::enable_if_t<
      containers::meta::is_array_pybind11<_rhs_type>::value
      > * = nullptr
    >
  OdeStorage(_rhs_type const & r)
    : auxRHS_{{_rhs_type(const_cast<_rhs_type &>(r).request())}}
  {}
#endif

  ~OdeStorage() = default;

  std::array<rhs_type, 1> auxRHS_;
};
//--------------------------------------------------

}}}//end namespace rompp::ode::impl
#endif
