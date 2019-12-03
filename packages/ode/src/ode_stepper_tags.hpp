/*
//@HEADER
// ************************************************************************
//
// ode_stepper_tags.hpp
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

#ifndef ODE_STEPPER_TAGS_HPP_
#define ODE_STEPPER_TAGS_HPP_

namespace pressio{ namespace ode{

namespace explicitmethods{
struct Undefined{};
struct Euler{};
struct RungeKutta4{};
}//end namespace pressio::ode::explicitmethods

namespace implicitmethods{
struct Undefined{};
struct Euler{};
struct BDF2{};

// this is used to define a stepper that is user-defined,
// so there is no specific info about it and it can be
// anything since the user assembles the time-discrete operators
struct Arbitrary{};

}//end namespace pressio::ode::implicitmethods


/* for now, support also enums for backward compability */
enum class ExplicitEnum
  {Undefined, Euler, RungeKutta4};

enum class ImplicitEnum
  {Undefined, Euler, BDF2,
   Arbitrary // this is used to define a stepper that is user-defined,
	     // so there is no specific info about it and it can be
             // anything since the user assembles the time-discrete operators
  };


namespace impl{

// helper to convert an enum to tag
template <::pressio::ode::ExplicitEnum>
struct ExplicitEnumToTagType{ using type = void; };

template <>
struct ExplicitEnumToTagType<::pressio::ode::ExplicitEnum::Undefined>{
  using type = ::pressio::ode::explicitmethods::Undefined;
};

template <>
struct ExplicitEnumToTagType<::pressio::ode::ExplicitEnum::Euler>{
  using type = ::pressio::ode::explicitmethods::Euler;
};

template <>
struct ExplicitEnumToTagType<::pressio::ode::ExplicitEnum::RungeKutta4>{
  using type = ::pressio::ode::explicitmethods::RungeKutta4;
};

// helper to convert an enum to tag
template <::pressio::ode::ImplicitEnum>
struct ImplicitEnumToTagType{ using type = void; };

template <>
struct ImplicitEnumToTagType<::pressio::ode::ImplicitEnum::Undefined>{
  using type = ::pressio::ode::implicitmethods::Undefined;
};

template <>
struct ImplicitEnumToTagType<::pressio::ode::ImplicitEnum::Euler>{
  using type = ::pressio::ode::implicitmethods::Euler;
};

template <>
struct ImplicitEnumToTagType<::pressio::ode::ImplicitEnum::BDF2>{
  using type = ::pressio::ode::implicitmethods::BDF2;
};

template <>
struct ImplicitEnumToTagType<::pressio::ode::ImplicitEnum::Arbitrary>{
  using type = ::pressio::ode::implicitmethods::Arbitrary;
};

}// end pressio::ode::impl

}}//end namespace pressio::ode
#endif
