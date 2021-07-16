/*
//@HEADER
// ************************************************************************
//
// ode_call_collector_dispatcher.hpp
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

#ifndef ODE_INTEGRATORS_IMPL_ODE_CALL_COLLECTOR_DISPATCHER_HPP_
#define ODE_INTEGRATORS_IMPL_ODE_CALL_COLLECTOR_DISPATCHER_HPP_

namespace pressio{ namespace ode{ namespace impl{

template <
  typename collector_type, typename time_type, typename state_type,
  typename enable = void
  >
struct CallCollectorDispatch;

// specialize for when collector accepts a native type
template <
  typename collector_type, typename time_type, typename state_type
  >
struct CallCollectorDispatch<
  collector_type, time_type, state_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_wrapper<state_type>::value and
    ::pressio::ode::constraints::collector_callable_with_step_time_native_container_return_void<
      collector_type, time_type, state_type
      >::value
    >
  >
{
  static void execute(collector_type	& collector,
		      const ::pressio::ode::types::step_t	& step,
		      const time_type	& time,
		      const state_type	& yIn){

    collector(step, time, *yIn.data());
  }
};

// specialize for when collector accepts a pressio container directly
template <
  typename collector_type, typename time_type, typename state_type
  >
struct CallCollectorDispatch<
  collector_type, time_type, state_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::predicates::is_wrapper<state_type>::value and
    ::pressio::ode::constraints::collector_callable_with_step_time_pressio_container_return_void<
      collector_type, time_type, state_type
      >::value
    >
  >
{
  static void execute(collector_type	& collector,
		      const ::pressio::ode::types::step_t	& step,
		      const time_type	& time,
		      const state_type	& yIn){

    collector(step, time, yIn);
  }
};

}}}//end namespace pressio::ode::impl
#endif  // ODE_INTEGRATORS_IMPL_ODE_CALL_COLLECTOR_DISPATCHER_HPP_
