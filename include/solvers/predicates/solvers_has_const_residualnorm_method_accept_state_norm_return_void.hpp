/*
//@HEADER
// ************************************************************************
//
// solvers_has_const_residualnorm_method_accept_state_norm_return_void.hpp
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

#ifndef SOLVERS_PREDICATES_SOLVERS_HAS_CONST_RESIDUALNORM_METHOD_ACCEPT_STATE_NORM_RETURN_VOID_HPP_
#define SOLVERS_PREDICATES_SOLVERS_HAS_CONST_RESIDUALNORM_METHOD_ACCEPT_STATE_NORM_RETURN_VOID_HPP_

namespace pressio{ namespace solvers{ namespace predicates {
  

template <
  typename T,
  typename state_t,
  typename norm_t,
  typename = void
  >
struct has_const_residualnorm_method_accept_state_norm_return_void
  : std::false_type{};

template <
  typename T,
  typename state_t,
  typename norm_t
  >
struct has_const_residualnorm_method_accept_state_norm_return_void<
  T, state_t, norm_t, 
  mpl::enable_if_t<
    std::is_void<
      decltype(
         std::declval<T const>().residualNorm
            (
              std::declval<state_t const &>(),
              ::pressio::Norm::Undefined,
              std::declval<norm_t &>()
            )
         )
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::solvers::predicates
#endif  // SOLVERS_PREDICATES_SOLVERS_HAS_CONST_RESIDUALNORM_METHOD_ACCEPT_STATE_NORM_RETURN_VOID_HPP_
