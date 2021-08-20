/*
//@HEADER
// ************************************************************************
//
// pressio_ode_common.hpp
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

#ifndef ODE_PRESSIO_ODE_COMMON_HPP_
#define ODE_PRESSIO_ODE_COMMON_HPP_

namespace pressio{ namespace ode{ 

struct ForwardEuler{};
struct RungeKutta4{};
struct AdamsBashforth2{};
struct SSPRungeKutta3{};

struct BDF1{};
struct BDF2{};
struct CrankNicolson{};
struct ImplicitArbitrary{};

template <typename T> struct is_explicit_stepper_tag : std::false_type{};
template <> struct is_explicit_stepper_tag<ForwardEuler> : std::true_type{};
template <> struct is_explicit_stepper_tag<RungeKutta4> : std::true_type{};
template <> struct is_explicit_stepper_tag<AdamsBashforth2> : std::true_type{};
template <> struct is_explicit_stepper_tag<SSPRungeKutta3> : std::true_type{};

template <typename T> struct is_implicit_stepper_tag : std::false_type{};
template <> struct is_implicit_stepper_tag<BDF1> : std::true_type{};
template <> struct is_implicit_stepper_tag<BDF2> : std::true_type{};
template <> struct is_implicit_stepper_tag<CrankNicolson> : std::true_type{};
template <> struct is_implicit_stepper_tag<ImplicitArbitrary> : std::true_type{};

// is_stepper_tag
template <typename T>
struct is_stepper_tag{
  static constexpr auto value = 
  is_explicit_stepper_tag<T>::value or is_implicit_stepper_tag<T>::value;
};

constexpr int implicit_stencil_size(::pressio::ode::BDF1){ return 2; }
constexpr int implicit_stencil_size(::pressio::ode::BDF2){ return 3; }
constexpr int implicit_stencil_size(::pressio::ode::CrankNicolson){ return 2; }

//! Default type for the order of a stepper
using stepper_order_type = int32_t;

// to set implicit stepper order for arbitrary stepper
template <stepper_order_type valueIn>
struct StepperOrder{
  static constexpr stepper_order_type value = valueIn;
};

// this is used to set the TOTAL number of states
// when the user chooses the arbitrary one
template <std::size_t valueIn>
struct StepperTotalNumberOfStates{ 
  static constexpr std::size_t value = valueIn; 
};

namespace constants{
template <typename scalar_t>
struct bdf1{
  using cnst = ::pressio::utils::Constants<scalar_t>;
  static constexpr scalar_t c_np1_= cnst::one();
  static constexpr scalar_t c_n_  = cnst::negOne();
  static constexpr scalar_t c_f_  = cnst::negOne();
};

template <typename scalar_t>
struct bdf2{
  using cnst = ::pressio::utils::Constants<scalar_t>;
  static constexpr scalar_t c_np1_ = cnst::one();
  static constexpr scalar_t c_n_   = cnst::negOne()*cnst::fourOvThree();
  static constexpr scalar_t c_nm1_ = cnst::oneOvThree();
  static constexpr scalar_t c_f_   = cnst::negOne()*cnst::twoOvThree();
};

template <typename scalar_t>
struct cranknicolson{
  using cnst = ::pressio::utils::Constants<scalar_t>;
  static constexpr scalar_t c_np1_  = cnst::one();
  static constexpr scalar_t c_n_    = cnst::negOne();
  static constexpr scalar_t c_fnp1_ = cnst::negOneHalf();
  static constexpr scalar_t c_fn_   = cnst::negOneHalf();
};
}//end namespace pressio::ode::constants

class nPlusOne{};
class n{};
class nMinusOne{};
class nMinusTwo{};
class nMinusThree{};
class nMinusFour{};

}}//end namespace pressio::ode

#include "./impl/ode_stencil_data_container.hpp"

namespace pressio{ namespace ode{ 

template<typename VelocityType, std::size_t N>
using ImplicitStencilVelocitiesContainer 
  = impl::StencilDataContainerImpl<VelocityType, N, nPlusOne /*stencil ends with n+1*/>;

template<typename StateType, std::size_t N>
using ImplicitStencilStatesContainer 
  = impl::StencilDataContainerImpl<StateType, N, n /*starts end at n*/>;

}}//end namespace pressio::ode

/*
   NOTE that the order below matters!
   Includes are ordered properly to avoid a tangled system.

   NOTE also that this header by itself means nothing and if you use
   it as such, you need to know what you are doing.
   This header is here to help the publicly-exposed includes named
   "pressio_ode_bla.hpp" inside the pressio/packages directory.
   Users of pressio should NOT rely on this file, but only
   on the top-level "pressio_ode_{explicit,implicit}.hpp".
*/
#include "ode_exceptions.hpp"
#include "predicates/ode_has_const_create_velocity_method_return_result.hpp"
#include "predicates/ode_has_const_velocity_method_accept_state_time_result_return_void.hpp"
#include "predicates/ode_has_const_create_discrete_time_residual_method_return_result.hpp"
#include "predicates/ode_has_const_discrete_time_residual_method_accept_step_time_dt_result_states_return_void.hpp"
#include "predicates/ode_has_const_create_discrete_time_jacobian_method_return_result.hpp"
#include "predicates/ode_has_const_discrete_time_jacobian_method_accepting_n_states_returning_void.hpp"
#include "predicates/ode_has_const_create_jacobian_method_return_result.hpp"
#include "predicates/ode_has_const_jacobian_method_accept_state_time_result_return_void.hpp"

#include "constraints/ode_continuous_time_system_with_at_least_velocity.hpp"
#include "constraints/ode_continuous_time_system_with_user_provided_jacobian.hpp"
#include "constraints/ode_discrete_time_system_with_user_provided_jacobian.hpp"
#include "constraints/ode_explicit_state.hpp"
#include "constraints/ode_explicit_velocity.hpp"
#include "constraints/ode_implicit_state.hpp"
#include "constraints/ode_implicit_residual.hpp"
#include "constraints/ode_implicit_jacobian.hpp"
#include "constraints/ode_implicit_residual_policy.hpp"
#include "constraints/ode_implicit_jacobian_policy.hpp"

#endif  // ODE_PRESSIO_ODE_COMMON_HPP_