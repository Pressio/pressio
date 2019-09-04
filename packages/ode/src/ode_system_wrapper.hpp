/*
//@HEADER
// ************************************************************************
//
// ode_system_wrapper.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef ODE_SYSTEM_WRAPPER_HPP_
#define ODE_SYSTEM_WRAPPER_HPP_

#include "ode_ConfigDefs.hpp"
#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#endif

namespace pressio{ namespace ode{ namespace impl{

template<typename model_type>
struct OdeSystemWrapper<
  model_type
#ifdef HAVE_PYBIND11
  , mpl::enable_if_t<
      ::pressio::mpl::not_same<model_type, pybind11::object >::value
      >
#endif
  >
{
  OdeSystemWrapper(const model_type & system)
    : data_(system){}

  OdeSystemWrapper() = delete;
  ~OdeSystemWrapper() = default;

  const model_type & get() const{
    return data_;
  }

private:
  const model_type & data_;
};


/* for some reason to be determined, when we deal with
 * python objects, we need to pass by copy
 */
#ifdef HAVE_PYBIND11
template<typename model_type>
struct OdeSystemWrapper<
  model_type,
  mpl::enable_if_t<
    ::pressio::mpl::is_same<model_type, pybind11::object >::value
    >
  >
{
  OdeSystemWrapper(const model_type & system)
    : data_(system){}

  OdeSystemWrapper() = delete;
  ~OdeSystemWrapper() = default;

  const model_type & get() const{
    return data_;
  }

private:
    model_type data_;
};
#endif

}}}//end namespace pressio::ode::impl
#endif
