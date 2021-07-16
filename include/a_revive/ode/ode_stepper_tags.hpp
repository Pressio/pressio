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


#ifndef ODE_ODE_STEPPER_TAGS_HPP_
#define ODE_ODE_STEPPER_TAGS_HPP_

namespace pressio{ namespace ode{

namespace explicitmethods{
struct Undefined{};
struct Euler{};
struct RungeKutta4{};
struct AdamsBashforth2{};
}//end namespace explicitmethods

namespace implicitmethods{
struct Undefined{};
struct BDF1{};
struct BDF2{};
struct CrankNicolson{};
using Euler = BDF1;

// this is used to define a stepper that is user-defined,
// so there is no specific info about it and it can be
// anything since the user assembles the time-discrete operators
struct Arbitrary{};
}//end namespace implicitmethods


// makes sense to leave the metafunctions here because every time
// a new tag is added, a corresponding metaf is needed too
namespace predicates{
// explicit
template <typename T>
struct is_explicit_stepper_tag : std::false_type{};

template <>
struct is_explicit_stepper_tag<explicitmethods::Euler> : std::true_type{};

template <>
struct is_explicit_stepper_tag<explicitmethods::RungeKutta4> : std::true_type{};

template <>
struct is_explicit_stepper_tag<explicitmethods::AdamsBashforth2> : std::true_type{};

// implicit
template <typename T>
struct is_implicit_stepper_tag : std::false_type{};

template <>
struct is_implicit_stepper_tag<implicitmethods::BDF1> : std::true_type{};

template <>
struct is_implicit_stepper_tag<implicitmethods::BDF2> : std::true_type{};

template <>
struct is_implicit_stepper_tag<implicitmethods::CrankNicolson> : std::true_type{};

template <>
struct is_implicit_stepper_tag<implicitmethods::Arbitrary> : std::true_type{};

// is_stepper_tag
template <typename T>
struct is_stepper_tag
{
  static constexpr auto value = is_explicit_stepper_tag<T>::value
    or is_implicit_stepper_tag<T>::value;
};
}//end namespace predicates

}}//end namespace pressio::ode
#endif  // ODE_ODE_STEPPER_TAGS_HPP_
