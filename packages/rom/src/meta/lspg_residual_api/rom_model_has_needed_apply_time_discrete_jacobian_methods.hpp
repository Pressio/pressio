/*
//@HEADER
// ************************************************************************
//
// rom_model_has_needed_apply_time_discrete_jacobian_methods.hpp
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

#ifndef ROM_MODEL_HAS_NEEDED_APPLY_TIME_DISCRETE_JACOBIAN_METHODS_HPP_
#define ROM_MODEL_HAS_NEEDED_APPLY_TIME_DISCRETE_JACOBIAN_METHODS_HPP_

#include "rom_has_apply_time_discrete_jacobian_method_accepting_n_states_returning_void.hpp"
#include "rom_has_apply_time_discrete_jacobian_method_accepting_n_states_returning_non_void.hpp"

namespace pressio{ namespace rom{ namespace meta {

template<
  typename model_t, typename step_t, typename scalar_t, typename state_t, typename dense_mat_t,
  typename enable = void
  >
struct model_has_needed_apply_time_discrete_jacobian_methods
  : std::false_type{};

template<
  typename model_t, typename step_t, typename scalar_t, typename state_t, typename dense_mat_t
  >
struct model_has_needed_apply_time_discrete_jacobian_methods<
  model_t, step_t, scalar_t, state_t, dense_mat_t,
  mpl::enable_if_t<
    has_apply_time_discrete_jacobian_method_accepting_n_states_returning_void<
      model_t, 2, step_t, scalar_t, state_t, dense_mat_t
      >::value and
    has_apply_time_discrete_jacobian_method_accepting_n_states_returning_void<
      model_t, 3, step_t, scalar_t, state_t, dense_mat_t
      >::value and
    has_apply_time_discrete_jacobian_method_accepting_n_states_returning_void<
      model_t, 4, step_t, scalar_t, state_t, dense_mat_t
      >::value and
    has_apply_time_discrete_jacobian_method_accepting_n_states_returning_non_void<
      model_t, 2, step_t, scalar_t, state_t, dense_mat_t
      >::value and
    has_apply_time_discrete_jacobian_method_accepting_n_states_returning_non_void<
      model_t, 3, step_t, scalar_t, state_t, dense_mat_t
      >::value and
    has_apply_time_discrete_jacobian_method_accepting_n_states_returning_non_void<
      model_t, 4, step_t, scalar_t, state_t, dense_mat_t
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::rom::meta
#endif
