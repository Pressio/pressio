/*
//@HEADER
// ************************************************************************
//
// ode_explicit_velocity_standard_policy.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef ODE_POLICIES_STANDARD_EXPLICIT_VELOCITY_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_EXPLICIT_VELOCITY_STANDARD_POLICY_HPP_

#include "../../ode_fwd.hpp"
#include "ode_explicit_velocity_policy_base.hpp"

#ifdef HAVE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#endif

namespace pressio{ namespace ode{ namespace policy{

/*
 * state_type = velocity_type
 * both are wrappers from containers
 */
template<
  typename state_type,
  typename system_type
  >
class ExplicitVelocityStandardPolicy<
  state_type, system_type, state_type,
  mpl::enable_if_t<
    containers::meta::is_vector_wrapper<state_type>::value
#ifdef HAVE_PYBIND11
    and mpl::not_same<system_type, pybind11::object >::value
#endif
    >
  >
  : public ExplicitVelocityPolicyBase<
  ExplicitVelocityStandardPolicy<state_type, system_type> >{

  using base_t = ExplicitVelocityPolicyBase<
    ExplicitVelocityStandardPolicy<state_type, system_type>>;
  friend base_t;

public:
  ExplicitVelocityStandardPolicy() = default;
  ~ExplicitVelocityStandardPolicy() = default;

  template < typename scalar_type >
  void operator()(const state_type & y,
		  state_type & f,
		  const system_type & model,
		  scalar_type t) const{
    model.velocity(*y.data(), t, *f.data());
  }

  template < typename scalar_type >
  state_type operator()(const state_type & y,
			const system_type & model,
			scalar_type t) const{
    return state_type(model.velocity(*y.data(), t));
  }
};//end class



/*
 * state_type = velocity_type
 * both are pybind11::array_t
 */
#ifdef HAVE_PYBIND11
template<
  typename state_type,
  typename system_type
  >
class ExplicitVelocityStandardPolicy<
  state_type, system_type, state_type,
  mpl::enable_if_t<
    mpl::is_same<system_type, pybind11::object >::value and
    containers::meta::is_cstyle_array_pybind11<state_type>::value
    >
  >
  : public ExplicitVelocityPolicyBase<
  ExplicitVelocityStandardPolicy<state_type, system_type> >{

  using base_t = ExplicitVelocityPolicyBase<
    ExplicitVelocityStandardPolicy<state_type, system_type>>;
  friend base_t;

public:
  ExplicitVelocityStandardPolicy() = default;
  ~ExplicitVelocityStandardPolicy() = default;

  template <typename scalar_type>
  void operator()(const state_type & y,
		  state_type & f,
		  const system_type & model,
		  scalar_type t) const{
    //printf("C++ f address: %p\n", f.data());
    model.attr("velocity")(y, t, f);
  }

  template <typename scalar_type>
  state_type operator()(const state_type & y,
  			const system_type & model,
  			scalar_type t) const{
    return model.attr("velocity")(y, t);
  }

};//end class
#endif


}}}//end namespace pressio::ode::policy
#endif
