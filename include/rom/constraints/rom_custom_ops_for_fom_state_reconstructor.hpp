/*
//@HEADER
// ************************************************************************
//
// rom_custom_ops_for_fom_state_reconstructor.hpp
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

#ifndef ROM_CONSTRAINTS_ROM_CUSTOM_OPS_FOR_FOM_STATE_RECONSTRUCTOR_HPP_
#define ROM_CONSTRAINTS_ROM_CUSTOM_OPS_FOR_FOM_STATE_RECONSTRUCTOR_HPP_

namespace pressio{ namespace rom{ namespace constraints {

template<
  typename T,
  typename fom_state_type,
  typename enable = void
  >
struct custom_ops_for_fom_state_reconstructor
  : std::false_type{};

template <
  typename T,
  typename fom_state_type
  >
struct custom_ops_for_fom_state_reconstructor<
  T, fom_state_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_wrapper<fom_state_type>::value
    and
    ::pressio::ops::predicates::has_method_deep_copy<
      T,
      typename ::pressio::containers::details::traits<fom_state_type>::wrapped_t,
      typename ::pressio::containers::details::traits<fom_state_type>::wrapped_t
      >::value
    and
    ::pressio::ops::predicates::has_method_set_zero<
      T,
      typename ::pressio::containers::details::traits<fom_state_type>::wrapped_t
      >::value
    and
    ::pressio::ops::predicates::has_method_axpy<
      T,
      typename ::pressio::containers::details::traits<fom_state_type>::wrapped_t,
      typename ::pressio::containers::details::traits<fom_state_type>::wrapped_t,
      typename ::pressio::containers::details::traits<fom_state_type>::scalar_t
      >::value
    >
  > : std::true_type{};

}}} // namespace pressio::rom::constraints
#endif  // ROM_CONSTRAINTS_ROM_CUSTOM_OPS_FOR_FOM_STATE_RECONSTRUCTOR_HPP_
