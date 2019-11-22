/*
//@HEADER
// ************************************************************************
//
// rom_mask_decorator_residual.hpp
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

#ifndef ROM_MASK_DECORATOR_RESIDUAL_HPP_
#define ROM_MASK_DECORATOR_RESIDUAL_HPP_

#include "../rom_fwd.hpp"
#include "../../../ode/src/implicit/policies/meta/ode_is_legitimate_implicit_residual_policy.hpp"

namespace pressio{ namespace rom{ namespace decorator{

template <typename maskable>
class Masked<maskable,
  //enable when the maskable type is a residual
  ::pressio::mpl::enable_if_t<maskable::isResidualPolicy_>
  > : public maskable
{
  using typename maskable::residual_t;
  using maskable::R_;
  mutable std::shared_ptr<residual_t> maskedR_ = {};

public:
  Masked() = delete;
  Masked(const maskable & obj) : maskable(obj){}

  template <typename ... Args>
  Masked(Args && ... args)
    : maskable(std::forward<Args>(args)...){}

  ~Masked() = default;

public:

  template <
    ::pressio::ode::ImplicitEnum odeMethod,
    int n,
    typename ode_state_t,
    typename ode_res_t,
    typename app_t,
    typename scalar_t
  >
  void operator()(const ode_state_t & odeY,
  		  ode_res_t & R,
		  const ::pressio::ode::AuxStatesContainer<false, ode_state_t,n> & prevStates,
  		  const app_t & app,
		  const scalar_t & t,
		  const scalar_t & dt,
		  const ::pressio::ode::types::step_t & step) const
  {
    maskable::template operator()<
      odeMethod, n>(odeY, R_, prevStates, app, t, dt, step);

    app.applyMask(*R_.data(), *R.data(), t);
  }


  template <
    ode::ImplicitEnum odeMethod,
    int n,
    typename ode_state_t,
    typename app_t,
    typename scalar_t
    >
  residual_t operator()(const ode_state_t & odeY,
			const ::pressio::ode::AuxStatesContainer<false, ode_state_t,n> & prevStates,
			const app_t & app,
			const scalar_t & t,
			const scalar_t & dt,
			const ::pressio::ode::types::step_t & step) const
  {
    auto R1 = maskable::template operator()
      <odeMethod, n>(odeY, prevStates, app, t, dt, step);

    if (!maskedR_){
      maskedR_ = std::make_shared
	<residual_t>( app.applyMask(*R1.data(), t) );
    }
    else
      app.applyMask(*R1.data(), *maskedR_->data(), t);

    return *maskedR_;
  }

};//end class

}}} //end namespace pressio::rom::decorator
#endif
