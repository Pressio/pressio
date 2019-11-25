/*
//@HEADER
// ************************************************************************
//
// rom_mask_decorator_jacobian.hpp
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

#ifndef ROM_MASK_DECORATOR_JACOBIAN_HPP_
#define ROM_MASK_DECORATOR_JACOBIAN_HPP_

#include "../rom_fwd.hpp"
#include "../../../ode/src/implicit/policies/meta/ode_is_legitimate_implicit_jacobian_policy.hpp"

namespace pressio{ namespace rom{ namespace decorator{


template <typename maskable>
class Masked<maskable,
  /* enable when decorating a jacobian policy */
  ::pressio::mpl::enable_if_t<maskable::isResidualPolicy_ == false>
  > : public maskable{

public:
  using typename maskable::apply_jac_return_t;
  mutable std::shared_ptr<apply_jac_return_t> maskedJJ_ = {};
  using maskable::JJ_;

public:
  Masked() = delete;
  Masked(const maskable & obj) : maskable(obj){}

  template <typename ... Args>
  Masked(Args && ... args)
    : maskable(std::forward<Args>(args)...){}

  ~Masked() = default;

public:

  template <
    typename stepper_tag,
    typename ode_state_t,
    typename app_t,
    typename scalar_t,
    typename ode_jac_t
  >
  void operator()(const ode_state_t & odeY,
  		  const app_t & app,
		  const scalar_t & t,
		  const scalar_t & dt,
		  const ::pressio::ode::types::step_t & step,
		  ode_jac_t & odeJJ) const
  {
    maskable::template operator()<stepper_tag>(odeY, app, t, dt, step, JJ_);
    app.applyMask(*JJ_.data(), *odeJJ.data(), t);
  }

  template <
    typename stepper_tag,
    typename ode_state_t,
    typename app_t,
    typename scalar_t
    >
  apply_jac_return_t operator()(const ode_state_t & odeY,
				const app_t & app,
				const scalar_t & t,
				const scalar_t & dt,
				const ::pressio::ode::types::step_t & step) const
  {
    maskable::template operator()<stepper_tag>(odeY, app, t, dt, step, JJ_);
    if (!maskedJJ_){
      maskedJJ_ = std::make_shared<apply_jac_return_t>( app.applyMask(*JJ_.data(), t) );
    }
    else
      app.applyMask(*JJ_.data(), *maskedJJ_->data(), t);

    return *maskedJJ_;
  }


};//end class

}}} //end namespace pressio::rom::decorator
#endif
