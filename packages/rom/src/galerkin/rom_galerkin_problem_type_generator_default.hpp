/*
//@HEADER
// ************************************************************************
//
// rom_galerkin_problem_type_generator_default.hpp
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

#ifndef ROM_GALERKIN_PROBLEM_TYPE_GENERATOR_DEFAULT_HPP_
#define ROM_GALERKIN_PROBLEM_TYPE_GENERATOR_DEFAULT_HPP_

namespace pressio{ namespace rom{ namespace galerkin{

template <
  typename stepper_tag,
  typename galerkin_state_type,
  typename ...Args
  >
struct DefaultProblemType
  : CommonTypes<galerkin_state_type, Args...>
{

  using base_t = CommonTypes<galerkin_state_type, Args...>;

  using typename base_t::fom_t;
  using typename base_t::scalar_t;
  using typename base_t::fom_native_state_t;
  using typename base_t::fom_state_t;
  using typename base_t::fom_velocity_t;
  using typename base_t::galerkin_state_t;
  using typename base_t::galerkin_residual_t;
  using typename base_t::decoder_t;
  using typename base_t::decoder_jac_t;
  using typename base_t::fom_state_reconstr_t;
  using typename base_t::fom_states_data;
  using typename base_t::ud_ops_t;

  // policy for evaluating the ode velocity
  using galerkin_residual_policy_t =
    ::pressio::rom::galerkin::DefaultExplicitVelocityPolicy<
    fom_states_data, fom_velocity_t, decoder_t, ud_ops_t>;

  // declare type of stepper object
  using galerkin_stepper_t = ::pressio::ode::explicitmethods::Stepper<
    stepper_tag, galerkin_state_type, fom_t,
    galerkin_residual_t, galerkin_residual_policy_t, scalar_t>;

};//end class

}}}//end  namespace pressio::rom::galerkin
#endif
