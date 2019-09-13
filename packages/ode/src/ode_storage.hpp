/*
//@HEADER
// ************************************************************************
//
// ode_storage.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
