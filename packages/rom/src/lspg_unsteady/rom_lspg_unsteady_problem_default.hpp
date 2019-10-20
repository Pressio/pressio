/*
//@HEADER
// ************************************************************************
//
// rom_lspg_unsteady_problem_default.hpp
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

#ifndef ROM_LSPG_UNSTEADY_PROBLEM_TYPE_GENERATOR_DEFAULT_HPP_
#define ROM_LSPG_UNSTEADY_PROBLEM_TYPE_GENERATOR_DEFAULT_HPP_

#include "./impl_velocity_api/rom_lspg_unsteady_problem_type_generator_default_velocity_api.hpp"
#include "./impl_residual_api/rom_lspg_unsteady_problem_type_generator_default_residual_api.hpp"

namespace pressio{ namespace rom{

namespace impl{

template <typename T, typename enable = void>
struct DefaultLSPGUnsteadyHelper{
  template <::pressio::ode::ImplicitEnum name, typename lspg_state_t, typename ...Args>
  using type = void;
};

template <typename T>
struct DefaultLSPGUnsteadyHelper<
  T,
  mpl::enable_if_t<
    ::pressio::rom::meta::model_meets_velocity_api_for_unsteady_lspg<T>::value
    >
  >
{
  template <::pressio::ode::ImplicitEnum name, typename lspg_state_t, typename ...Args>
  using type = impl::DefaultLSPGUnsteadyTypeGeneratorVelocityAPI<name, T, lspg_state_t, Args...>;
};


template <typename T>
struct DefaultLSPGUnsteadyHelper<
  T,
  mpl::enable_if_t<
    ::pressio::rom::meta::model_meets_residual_api_for_unsteady_lspg<T>::value
    >
  >
{
  template <::pressio::ode::ImplicitEnum name, typename lspg_state_t, typename ...Args>
  using type = impl::DefaultLSPGUnsteadyTypeGeneratorResidualAPI<name, T, lspg_state_t, Args...>;
};

}// end namespace pressio::rom::impl

template <
  ::pressio::ode::ImplicitEnum name,
  typename fom_type,
  typename lspg_state_type,
  typename ...Args
  >
using DefaultLSPGUnsteady =
  typename impl::DefaultLSPGUnsteadyHelper<fom_type>::template type<name, lspg_state_type, Args...>;

}}//end  namespace pressio::rom
#endif
