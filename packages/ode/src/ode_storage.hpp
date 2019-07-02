
#ifndef ODE_STORAGE_HPP_
#define ODE_STORAGE_HPP_

#include "ode_ConfigDefs.hpp"
#include <array>
#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#endif

namespace pressio{ namespace ode{ namespace impl{

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

template<typename T, int n>
struct OdeStorage;


/* n = 1 */
template<typename T>
struct OdeStorage<T, 1>{

  template <
    typename _T = T
#ifdef HAVE_PYBIND11
    , mpl::enable_if_t<
      !containers::meta::is_array_pybind11<_T>::value
      > * = nullptr
#endif
    >
  OdeStorage(_T const & y)
    : data_{{y}}{}

#ifdef HAVE_PYBIND11
  template <
    typename _T = T,
    mpl::enable_if_t<
      containers::meta::is_array_pybind11<_T>::value
      > * = nullptr
    >
  OdeStorage(_T const & y)
    : data_{{_T(const_cast<_T &>(y).request())}}
  {}
#endif

  OdeStorage() = delete;
  ~OdeStorage() = default;

  std::array<T, 1> data_;
};


/* n = 2 */
template<typename T>
struct OdeStorage<T, 2>{

  template <
    typename _T = T
#ifdef HAVE_PYBIND11
    , mpl::enable_if_t<
      !containers::meta::is_array_pybind11<_T>::value
      > * = nullptr
#endif
    >
  OdeStorage(_T const & y)
    : data_{{y,y}}{}

#ifdef HAVE_PYBIND11
  template <
    typename _T = T,
    mpl::enable_if_t<
      containers::meta::is_array_pybind11<_T>::value
      > * = nullptr
    >
  OdeStorage(_T const & y)
    : data_{{_T(const_cast<_T &>(y).request()),
	     _T(const_cast<_T &>(y).request())}}
  {}
#endif

  OdeStorage() = delete;
  ~OdeStorage() = default;

  std::array<T, 2> data_;
};


/* n = 3 */
template<typename T>
struct OdeStorage<T, 3>{

  template <
    typename _T = T
#ifdef HAVE_PYBIND11
    , mpl::enable_if_t<
      !containers::meta::is_array_pybind11<_T>::value
      > * = nullptr
#endif
    >
  OdeStorage(_T const & y)
    : data_{{y,y,y}}{}

#ifdef HAVE_PYBIND11
  template <
    typename _T = T,
    mpl::enable_if_t<
      containers::meta::is_array_pybind11<_T>::value
      > * = nullptr
    >
  OdeStorage(_T const & y)
    : data_{{_T(const_cast<_T &>(y).request()),
	     _T(const_cast<_T &>(y).request()),
	     _T(const_cast<_T &>(y).request())}}
  {}
#endif

  OdeStorage() = delete;
  ~OdeStorage() = default;

  std::array<T, 3> data_;
};


/* n = 4 */
template<typename T>
struct OdeStorage<T, 4>{

  template <
    typename _T = T
#ifdef HAVE_PYBIND11
    , mpl::enable_if_t<
      !containers::meta::is_array_pybind11<_T>::value
      > * = nullptr
#endif
    >
  OdeStorage(_T const & y)
    : data_{{y,y,y,y}}{}

#ifdef HAVE_PYBIND11
  template <
    typename _T = T,
    mpl::enable_if_t<
      containers::meta::is_array_pybind11<_T>::value
      > * = nullptr
    >
  OdeStorage(_T const & y)
    : data_{{_T(const_cast<_T &>(y).request()),
	     _T(const_cast<_T &>(y).request()),
	     _T(const_cast<_T &>(y).request()),
	     _T(const_cast<_T &>(y).request())}}
  {}
#endif

  OdeStorage() = delete;
  ~OdeStorage() = default;

  std::array<T, 4> data_;
};

}}}//end namespace pressio::ode::impl
#endif
