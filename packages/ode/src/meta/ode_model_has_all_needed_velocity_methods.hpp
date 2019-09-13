/*
//@HEADER
// ************************************************************************
//
// ode_model_has_all_needed_velocity_methods.hpp
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

#ifndef ODE_MODEL_HAS_ALL_NEEDED_VELOCITY_METHODS_HPP_
#define ODE_MODEL_HAS_ALL_NEEDED_VELOCITY_METHODS_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../../../mpl/src/detection_idiom.hpp"

namespace pressio{ namespace ode{ namespace meta {


template <typename T, typename state_t, typename sc_t, typename = void>
struct has_velocity_method_callable_with_two_args : std::false_type{};

template <typename T, typename state_t, typename sc_t>
struct has_velocity_method_callable_with_two_args<
  T, state_t, sc_t,
  ::pressio::mpl::enable_if_t<
    !std::is_void<
      decltype(
	       std::declval<T>().velocity(
					  std::declval<state_t const&>(),
					  std::declval<sc_t>()
					  )
	       )
      >::value
    >
  > : std::true_type{};



template <typename T,
	  typename state_t,
	  typename sc_t,
	  typename velo_t,
	  typename = void>
struct has_velocity_method_callable_with_three_args : std::false_type{};

template <typename T,
	  typename state_t,
	  typename sc_t,
	  typename velo_t>
struct has_velocity_method_callable_with_three_args<
  T, state_t, sc_t, velo_t,
  ::pressio::mpl::enable_if_t<
    std::is_void<
      decltype(
	       std::declval<T>().velocity(
					  std::declval<state_t const&>(),
					  std::declval<sc_t>(),
					  std::declval<velo_t &>()
					  )
	   )
      >::value
    >
  > : std::true_type{};

//---------------------------------------------------------------


template<
  typename model_type,
  typename state_type,
  typename velocity_type,
  typename scalar_type,
  typename enable = void
  >
struct model_has_needed_velocity_methods : std::false_type{};

template<
  typename model_type,
  typename state_type,
  typename velocity_type,
  typename scalar_type
  >
struct model_has_needed_velocity_methods<
  model_type, state_type, velocity_type, scalar_type,
  mpl::enable_if_t<
    // has method with 2 arguments,
    has_velocity_method_callable_with_two_args<
      model_type, state_type, scalar_type
      >::value and
    // has velocity method with 3 arguments
    has_velocity_method_callable_with_three_args<
      model_type, state_type, scalar_type, velocity_type
      >::value and
    // method with 2 arguments returns a velocity_type
    mpl::is_same<
      velocity_type,
      decltype(
	       std::declval<model_type>().velocity
	       (
		std::declval<state_type const&>(),
		std::declval<scalar_type>())
	       )
      >::value
    >
  > : std::true_type{};


}}} // namespace pressio::ode::meta
#endif
