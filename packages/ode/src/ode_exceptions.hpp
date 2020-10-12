/*
//@HEADER
// ************************************************************************
//
// ode_exceptions.hpp
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

#ifndef ODE_ODE_EXCEPTIONS_HPP_
#define ODE_ODE_EXCEPTIONS_HPP_

#include <exception>

namespace pressio{ namespace eh{

class time_step_failure
  : public std::exception
{
  std::string myerr_ = "Time step failed";
  std::string append_ = {};

public:
  time_step_failure() = default;
  time_step_failure(const time_step_failure &) = default;
  time_step_failure & operator=(const time_step_failure &) = default;
  time_step_failure(time_step_failure &&) = default;
  time_step_failure & operator=(time_step_failure &&) = default;
  ~time_step_failure() = default;

  explicit time_step_failure(std::string append)
    : append_{append}{
    myerr_ += append_;
  }

  const char * what () const throw (){
    return myerr_.c_str();
   }
};


class velocity_failure_unrecoverable
  : public std::exception
{
  std::string myerr_ = "Velocity evaluation failed";
  std::string append_ = {};

public:
  velocity_failure_unrecoverable() = default;
  velocity_failure_unrecoverable(const velocity_failure_unrecoverable &) = default;
  velocity_failure_unrecoverable & operator=(const velocity_failure_unrecoverable &) = default;
  velocity_failure_unrecoverable(velocity_failure_unrecoverable &&) = default;
  velocity_failure_unrecoverable & operator=(velocity_failure_unrecoverable &&) = default;
  ~velocity_failure_unrecoverable() = default;

  explicit velocity_failure_unrecoverable(std::string append)
    : append_{append}{
    myerr_ += append_;
  }

  const char * what () const throw (){
    return myerr_.c_str();
   }
};


class discrete_time_residual_failure_unrecoverable
  : public std::exception
{
  std::string myerr_ = "discreteTimeResidual failed";
  std::string append_ = {};

public:
  discrete_time_residual_failure_unrecoverable() = default;
  discrete_time_residual_failure_unrecoverable(const discrete_time_residual_failure_unrecoverable &) = default;
  discrete_time_residual_failure_unrecoverable & operator=(const discrete_time_residual_failure_unrecoverable &) = default;
  discrete_time_residual_failure_unrecoverable(discrete_time_residual_failure_unrecoverable &&) = default;
  discrete_time_residual_failure_unrecoverable & operator=(discrete_time_residual_failure_unrecoverable &&) = default;
  ~discrete_time_residual_failure_unrecoverable() = default;

  explicit discrete_time_residual_failure_unrecoverable(std::string append)
    : append_{append}{
    myerr_ += append_;
  }

  const char * what () const throw (){
    return myerr_.c_str();
   }
};

}}//end namespace pressio::eh
#endif  // ODE_ODE_EXCEPTIONS_HPP_
