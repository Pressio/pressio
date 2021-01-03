/*
//@HEADER
// ************************************************************************
//
// ode_required_number_of_states.hpp
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

#ifndef ODE_ODE_REQUIRED_NUMBER_OF_STATES_HPP_
#define ODE_ODE_REQUIRED_NUMBER_OF_STATES_HPP_

namespace pressio{ namespace ode{

// base case
template <typename stepper_tag>
struct requiredNumberOfStates
{
  static constexpr std::size_t value = 1;
};

template <>
struct requiredNumberOfStates<::pressio::ode::implicitmethods::Euler>
{
  // need to store: state_n+1, state_n
  static constexpr std::size_t value = 2;
};

template <>
struct requiredNumberOfStates<::pressio::ode::implicitmethods::BDF2>
{
  // need: state_n+1, state_n, state_n-1
  static constexpr std::size_t value = 3;
};

template <>
struct requiredNumberOfStates<::pressio::ode::implicitmethods::CrankNicolson>
{
  // need: state_n+1, state_n
  static constexpr std::size_t value = 2;
};

template <>
struct requiredNumberOfStates<::pressio::ode::explicitmethods::Euler>
{
  static constexpr std::size_t value = 1;
};

template <>
struct requiredNumberOfStates<::pressio::ode::explicitmethods::RungeKutta4>
{
  static constexpr std::size_t value = 1;
};

}} //end namespace pressio::ode
#endif  // ODE_ODE_REQUIRED_NUMBER_OF_STATES_HPP_
