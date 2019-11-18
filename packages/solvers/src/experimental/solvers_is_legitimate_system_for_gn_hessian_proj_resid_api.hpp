/*
//@HEADER
// ************************************************************************
//
// solvers_is_legitimate_system_for_gn_hessian_proj_resid_api.hpp
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

#ifndef SOLVERS_IS_LEGITIMATE_SYSTEM_FOR_GN_HESSIAN_PROJ_RESID_API_HPP_
#define SOLVERS_IS_LEGITIMATE_SYSTEM_FOR_GN_HESSIAN_PROJ_RESID_API_HPP_

#include "../meta/solvers_basic_meta.hpp"
#include "solvers_system_has_needed_compute_hessian_and_projected_residual.hpp"

namespace pressio{ namespace solvers{ namespace meta { namespace experimental{

template<typename system_type, typename enable = void>
struct is_legitimate_system_for_gn_hessian_proj_residual_api : std::false_type{};

template<typename system_type>
struct is_legitimate_system_for_gn_hessian_proj_residual_api
<
  system_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::mpl::is_detected<has_scalar_typedef, system_type>::value   and
    ::pressio::mpl::is_detected<has_state_typedef, system_type>::value    and
    /* check that system contains the hessian and j^T R, however that is called*/
    //
    system_has_needed_compute_hessian_and_proejcted_residual_methods</*...*/>::value
    >
  > : std::true_type{};

}}}} // namespace pressio::solvers::meta::experimental
#endif
