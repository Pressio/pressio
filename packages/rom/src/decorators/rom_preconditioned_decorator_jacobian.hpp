/*
//@HEADER
// ************************************************************************
//
// rom_preconditioned_decorator_jacobian.hpp
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

#ifndef ROM_PRECONDITIONED_DECORATOR_JACOBIAN_HPP_
#define ROM_PRECONDITIONED_DECORATOR_JACOBIAN_HPP_

#include "../rom_fwd.hpp"
#include "../../../ode/src/implicit/policies/meta/ode_is_legitimate_implicit_jacobian_policy.hpp"

namespace pressio{ namespace rom{ namespace decorator{

/* overload when decorating a jacobian policy */
template <typename preconditionable>
class Preconditioned<
  preconditionable,
  ::pressio::mpl::enable_if_t<preconditionable::isResidualPolicy_ == false>
  > : public preconditionable{

public:
  using typename preconditionable::apply_jac_return_t;
  using preconditionable::fomStates_;

public:
  Preconditioned() = delete;
  Preconditioned(const preconditionable & obj) : preconditionable(obj){}

  template <typename ... Args>
  Preconditioned(Args && ... args)
    : preconditionable(std::forward<Args>(args)...){}

  ~Preconditioned() = default;

public:

  //-------------------------------
  // for unsteady LSPG
  //-------------------------------
  template <
    ode::ImplicitEnum odeMethod,
    typename ode_state_t,
    typename app_t,
    typename scalar_t
  >
  apply_jac_return_t operator()(const ode_state_t & odeY,
				const app_t & app,
				scalar_t t, scalar_t dt) const
  {
    auto JJ = preconditionable::template operator()<odeMethod>(odeY, app, t, dt);
    const auto & yFom = fomStates_.getCRefToFomState();
    app.applyPreconditioner(*yFom.data(), *JJ.data(), t);
    return JJ;
  }

  template <
      ode::ImplicitEnum odeMethod,
      typename ode_state_t,
      typename ode_jac_t,
      typename app_t,
      typename scalar_t>
  void operator()(const ode_state_t & odeY, ode_jac_t & odeJJ,
  		  const app_t & app, scalar_t t, scalar_t dt) const
  {
    preconditionable::template operator()<odeMethod>(odeY, odeJJ, app, t, dt);
    const auto & yFom = fomStates_.getCRefToFomState();
    app.applyPreconditioner(*yFom.data(), *odeJJ.data(), t);
  }


  //-------------------------------
  // for steady LSPG
  //-------------------------------
  template <
      typename lspg_state_t,
      typename app_t>
  apply_jac_return_t operator()(const lspg_state_t & romY,
                                const app_t & app) const
  {
    auto JJ = preconditionable::template operator()
	<lspg_state_t, app_t>(romY, app);

    const auto & yFom = fomStates_.getCRefToFomState();
    app.applyPreconditioner(*yFom.data(), *JJ.data());
    return JJ;
  }

  template <
      typename lspg_state_t,
      typename lspg_jac_t,
      typename app_t>
  void operator()(const lspg_state_t & romY,
                  lspg_jac_t & romJJ,
                  const app_t & app) const
  {
    preconditionable::template operator()
	  <lspg_state_t, lspg_jac_t, app_t>(romY, romJJ, app);

    const auto & yFom = fomStates_.getCRefToFomState();
    app.applyPreconditioner(*yFom.data(), *romJJ.data());
  }

};//end class

}}} //end namespace pressio::rom::decorator

#endif
