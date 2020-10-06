/*
//@HEADER
// ************************************************************************
//
// ode_system_wrapper.hpp
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

#ifndef ODE_ODE_SYSTEM_WRAPPER_HPP_
#define ODE_ODE_SYSTEM_WRAPPER_HPP_

namespace pressio{ namespace ode{ namespace impl{

template<
  typename system_type,
  typename enable = void
  >
struct OdeSystemWrapper;


template<typename system_type>
struct OdeSystemWrapper<
  system_type
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
  , mpl::enable_if_t<
      ::pressio::mpl::not_same<system_type, pybind11::object >::value
      >
#endif
  >
{
  OdeSystemWrapper() = delete;

  OdeSystemWrapper(const system_type & system)
    : data_(system){}

  OdeSystemWrapper(const OdeSystemWrapper &) = default;
  OdeSystemWrapper & operator=(const OdeSystemWrapper &) = default;

  OdeSystemWrapper(OdeSystemWrapper &&) = default;
  OdeSystemWrapper & operator=(OdeSystemWrapper &&) = default;

  ~OdeSystemWrapper() = default;

  const system_type & get() const{
    return data_;
  }

private:
  std::reference_wrapper<const system_type> data_;
};


/* for some reason to be determined, when we deal with
 * python objects, we need to pass by copy
 */
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template<typename system_type>
struct OdeSystemWrapper<
  system_type,
  mpl::enable_if_t<
    ::pressio::mpl::is_same<system_type, pybind11::object >::value
    >
  >
{
  OdeSystemWrapper() = delete;

  OdeSystemWrapper(const system_type & system) : data_(system){}

  OdeSystemWrapper(const OdeSystemWrapper &) = default;
  OdeSystemWrapper & operator=(const OdeSystemWrapper &) = default;

  OdeSystemWrapper(OdeSystemWrapper &&) = default;
  OdeSystemWrapper & operator=(OdeSystemWrapper &&) = default;

  ~OdeSystemWrapper() = default;

  const system_type & get() const{
    return data_;
  }

private:
    system_type data_;
};
#endif

}}}//end namespace pressio::ode::impl
#endif  // ODE_ODE_SYSTEM_WRAPPER_HPP_
