/*
//@HEADER
// ************************************************************************
//
// exceptions.hpp
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

#ifndef PRESSIO_ODE_ODE_STRONG_TYPES_HPP_
#define PRESSIO_ODE_ODE_STRONG_TYPES_HPP_

namespace pressio{ namespace ode{

template<class T>
struct StepStartAt{
  using value_type = T;
  value_type value_{};
  StepStartAt() = default;
  constexpr explicit StepStartAt(value_type valueIn) : value_(std::move(valueIn)){}
  constexpr value_type get() const { return value_; }
};

template<class T>
struct StepEndAt{
  using value_type = T;
  value_type value_{};
  StepEndAt() = default;
  constexpr explicit StepEndAt(value_type valueIn) : value_(std::move(valueIn)){}
  constexpr value_type get() const { return value_; }
};

template<class T>
struct StepSize{
  using value_type = T;
  value_type value_{};
  StepSize() = default;
  constexpr explicit StepSize(value_type valueIn) : value_(std::move(valueIn)){}
  constexpr value_type get() const { return value_; }

  StepSize& operator= (const value_type& newVal){
    value_ = newVal;
    return *this;
  }
};

template<class T>
struct StepSizeMinAllowedValue{
  using value_type = T;
  value_type value_{};
  StepSizeMinAllowedValue() = default;
  constexpr explicit StepSizeMinAllowedValue(value_type valueIn) : value_(std::move(valueIn)){}
  constexpr value_type get() const { return value_; }

  StepSizeMinAllowedValue& operator= (const value_type& newVal){
    value_ = newVal;
    return *this;
  }
};

template<class T>
struct StepSizeScalingFactor{
  using value_type = T;
  value_type value_{};
  StepSizeScalingFactor() = default;
  constexpr explicit StepSizeScalingFactor(value_type valueIn) : value_(std::move(valueIn)){}
  constexpr value_type get() const { return value_; }

  StepSizeScalingFactor& operator= (const value_type& newVal){
    value_ = newVal;
    return *this;
  }
};

struct StepCount{
  using value_type = int32_t;
  value_type value_{};
  StepCount() = default;
  constexpr explicit StepCount(value_type valueIn) : value_(valueIn){}
  constexpr value_type get() const { return value_; }
};

struct IntermediateStepCount{
  using value_type = int32_t;
  value_type value_{};
  IntermediateStepCount() = default;
  constexpr explicit IntermediateStepCount(value_type valueIn) : value_(valueIn){}
  constexpr value_type get() const { return value_; }
};

}}// end namespace ode

#endif  // PRESSIO_ODE_ODE_STRONG_TYPES_HPP_
