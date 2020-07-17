/*
//@HEADER
// ************************************************************************
//
// rom_wls_hessian_and_gradient_sequential_policy_helper_impl.hpp
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

#ifndef ROM_WLS_IMPL_POLICIES_ROM_WLS_HESSIAN_AND_GRADIENT_SEQUENTIAL_POLICY_HELPER_IMPL_HPP_
#define ROM_WLS_IMPL_POLICIES_ROM_WLS_HESSIAN_AND_GRADIENT_SEQUENTIAL_POLICY_HELPER_IMPL_HPP_

#include "rom_wls_hessian_and_gradient_sequential_policy_impl.hpp"

namespace pressio{ namespace rom{ namespace wls{ namespace impl{

template<
  typename fom_system_type,
  typename decoder_type,
  typename ode_tag,
  typename hess_structure_tag,
  typename ...Args
  >
struct HessGradSeqPolHelper
{
  // find if Args contain a valid updating tag for jacobian, if not default = NonFrozen
  using jacUpdatingTagIC  = ::pressio::mpl::variadic::find_if_unary_pred_t<
    ::pressio::rom::wls::predicates::is_legitimate_jacobian_updating_tag, Args...>;
  using jac_update_t = ::pressio::mpl::variadic::at_or_t<::pressio::rom::wls::NonFrozenJacobian,
							  jacUpdatingTagIC::value, Args...>;

  // set the container type based on the tag
  using jac_t	   = typename decoder_type::jacobian_type;
  using jac_cont_t = typename std::conditional<
    std::is_same<jac_update_t, ::pressio::rom::wls::NonFrozenJacobian>::value,
    ::pressio::rom::wls::NonFrozenJacobiansContainer<jac_t>,
    ::pressio::rom::wls::FrozenJacobiansContainer<jac_t>
    >::type;

  // find if Args contain a valid precond type, if not default = NoPreconditioner
  using precIC  = ::pressio::mpl::variadic::find_if_unary_pred_t<::pressio::rom::wls::predicates::is_legitimate_preconditioner_type, Args...>;
  using prec_t = ::pressio::mpl::variadic::at_or_t<::pressio::rom::wls::preconditioners::NoPreconditioner, precIC::value, Args...>;

  // final type
  using type = ::pressio::rom::wls::impl::HessianGradientSequentialPolicy<
    fom_system_type, decoder_type, ode_tag, hess_structure_tag, prec_t, jac_cont_t>;
};

}}}} //end namespace pressio::rom::wls::impl
#endif  // ROM_WLS_IMPL_POLICIES_ROM_WLS_HESSIAN_AND_GRADIENT_SEQUENTIAL_POLICY_HELPER_IMPL_HPP_
