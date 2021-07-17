/*
//@HEADER
// ************************************************************************
//
// rom_compose_impl.hpp
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

#ifndef ROM_GALERKIN_IMPL_DISCRETE_TIME_API_ROM_COMPOSE_IMPL_HPP_
#define ROM_GALERKIN_IMPL_DISCRETE_TIME_API_ROM_COMPOSE_IMPL_HPP_

#include "../rom_problem_tags.hpp"
#include "../rom_galerkin_types_selector.hpp"

#include "./policies/rom_galerkin_fom_residual_policy.hpp"
#include "./policies/rom_galerkin_fom_apply_jacobian_policy.hpp"

#include "./traits/rom_galerkin_common_traits.hpp"
#include "./traits/rom_galerkin_default_problem_traits.hpp"
#include "./traits/rom_galerkin_masked_residual_problem_traits.hpp"
#include "./traits/rom_galerkin_hyper_reduced_residual_problem_traits.hpp"
#include "./rom_galerkin_default_problem.hpp"
#include "./rom_galerkin_masked_residual_problem.hpp"
#include "./rom_galerkin_hyper_reduced_residual_problem.hpp"

namespace pressio{ namespace rom{ namespace galerkin{ namespace impl{ namespace disctime{

template<
  typename problem_tag,
  typename dummy,
  typename stepper_tag,
  typename fom_system_t,
  typename decoder_type,
  typename galerkin_state_type,
  typename galerkin_jacobian_type,
  typename ...Args>
struct compose
{
  static constexpr bool api_error =
    ::pressio::rom::why_not_discrete_time_system_with_user_provided_apply_jacobian<
    fom_system_t, typename decoder_type::jacobian_type>::value;
  static_assert(api_error, "");

  using type = void;
};

/***
    default Galerkin
 ***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename galerkin_state_type,
  typename galerkin_jacobian_type,
  typename ...Args
  >
struct compose<
  ::pressio::rom::galerkin::impl::Default,
  mpl::enable_if_t<
    ::pressio::rom::constraints::discrete_time_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  stepper_tag, fom_system_type, decoder_type, galerkin_state_type, galerkin_jacobian_type, Args...>
{
  static_assert
  (std::is_same< stepper_tag, ::pressio::ode::implicitmethods::Arbitrary>::value,
   "Galerkin with discrete-time API currently only accepts Arbitrary stepper");

  using type = ::pressio::rom::galerkin::impl::DefaultProblemDiscreteTimeApi<
    stepper_tag, fom_system_type, galerkin_state_type,
    galerkin_jacobian_type, decoder_type, Args...>;
};

/***
    masked residual Galerkin
 ***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename galerkin_state_type,
  typename galerkin_jacobian_type,
  typename masker_type,
  typename projector_type,
  typename ...Args
  >
struct compose<
  ::pressio::rom::galerkin::impl::MaskedResidual,
  mpl::enable_if_t<
    ::pressio::rom::constraints::discrete_time_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  stepper_tag, fom_system_type, decoder_type, galerkin_state_type, galerkin_jacobian_type,
  masker_type, projector_type, Args...>
{
  static_assert
  (std::is_same< stepper_tag, ::pressio::ode::implicitmethods::Arbitrary>::value,
   "Galerkin with discrete-time API currently only accepts Arbitrary stepper");

  using type = ::pressio::rom::galerkin::impl::MaskedResidualProblemDiscreteTimeApi<
    stepper_tag, fom_system_type, galerkin_state_type, galerkin_jacobian_type,
    decoder_type, masker_type, projector_type, Args...>;
};

/***
    hyper-reduced residual Galerkin
 ***/
template<
  typename stepper_tag,
  typename fom_system_type,
  typename decoder_type,
  typename galerkin_state_type,
  typename galerkin_jacobian_type,
  typename projector_type,
  typename ...Args
  >
struct compose<
  ::pressio::rom::galerkin::impl::HyperReducedResidual,
  mpl::enable_if_t<
    ::pressio::rom::constraints::discrete_time_system_with_user_provided_apply_jacobian<
      fom_system_type, typename decoder_type::jacobian_type>::value
    >,
  stepper_tag, fom_system_type, decoder_type, galerkin_state_type, galerkin_jacobian_type,
  projector_type, Args...>
{
  static_assert
  (std::is_same< stepper_tag, ::pressio::ode::implicitmethods::Arbitrary>::value,
   "Galerkin with discrete-time API currently only accepts Arbitrary stepper");

  using type = ::pressio::rom::galerkin::impl::HyperReducedResidualProblemDiscreteTimeApi<
    stepper_tag, fom_system_type, galerkin_state_type, galerkin_jacobian_type,
    decoder_type, projector_type, Args...>;
};

//-------------------------------------------------------
} // end namespace pressio::rom::galerkin::impl::disctime
//-------------------------------------------------------

// default
template<typename ...Args>
using composeDefaultProblemDiscTime =
  impl::disctime::compose<
  impl::Default, void,
  typename std::remove_cv<typename std::remove_reference<Args>::type>::type...
  >;

template<typename ...Args>
using composeDefaultProblemDiscTime_t =
  typename composeDefaultProblemDiscTime<Args...>::type;

// masked residual
template<typename ...Args>
using composeMaskedResidualProblemDiscTime =
  impl::disctime::compose<
  impl::MaskedResidual, void,
  typename std::remove_cv<typename std::remove_reference<Args>::type>::type...
  >;

template<typename ...Args>
using composeMaskedResidualProblemDiscTime_t =
  typename composeMaskedResidualProblemDiscTime<Args...>::type;

// hyperReduced residual
template<typename ...Args>
using composeHyperReducedResidualProblemDiscTime =
  impl::disctime::compose<
  impl::HyperReducedResidual, void,
  typename std::remove_cv<typename std::remove_reference<Args>::type>::type...
  >;

template<typename ...Args>
using composeHyperReducedResidualProblemDiscTime_t =
  typename composeHyperReducedResidualProblemDiscTime<Args...>::type;

}}}} // end namespace pressio::rom::galerkin::impl
#endif  // ROM_GALERKIN_IMPL_DISCRETE_TIME_API_ROM_COMPOSE_IMPL_HPP_
