/*
//@HEADER
// ************************************************************************
//
// rom_preconditioned_decorator.hpp
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

#ifndef ROM_PRECONDITIONED_DECORATOR_HPP_
#define ROM_PRECONDITIONED_DECORATOR_HPP_

#include "../rom_fwd.hpp"
#include "../../../ode/src/implicit/policies/meta/ode_is_legitimate_implicit_residual_policy.hpp"
#include "../../../ode/src/implicit/policies/meta/ode_is_legitimate_implicit_jacobian_policy.hpp"

namespace pressio{ namespace rom{ namespace decorator{

/* overload when decorating a residual policy */
template <typename preconditionable>
class Preconditioned<
  preconditionable,
  ::pressio::mpl::enable_if_t<preconditionable::isResidualPolicy_>
  > : public preconditionable{

  using typename preconditionable::fom_rhs_t;
  using preconditionable::yFom_;

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
    int n,
    typename ode_state_t,
    typename app_t,
    typename scalar_t
  >
  fom_rhs_t operator()(const ode_state_t & odeY,
			 const std::array<ode_state_t, n> & oldYs,
			 const app_t & app,
			 scalar_t t,
			 scalar_t dt) const
  {
    auto result = preconditionable::template operator()<
      odeMethod, n>(odeY, oldYs, app, t, dt);

    app.applyPreconditioner(*yFom_.data(), *result.data(), t);
    return result;
  }

  template <
    ode::ImplicitEnum odeMethod,
    int n,
    typename ode_state_t,
    typename ode_res_t,
    typename app_t,
    typename scalar_t
    >
  void operator()(const ode_state_t & odeY,
  		  ode_res_t & odeR,
  		  const std::array<ode_state_t, n> & oldYs,
  		  const app_t & app,
		  scalar_t t,
		  scalar_t dt) const
  {
    preconditionable::template operator()<
      odeMethod, n>(odeY, odeR, oldYs, app, t, dt);

    app.applyPreconditioner(*yFom_.data(), *odeR.data(), t);
  }


  //-------------------------------
  // for steady LSPG
  //-------------------------------
  template <
      typename lspg_state_t,
      typename fom_t>
  fom_rhs_t operator()(const lspg_state_t  & romY,
                         const fom_t   & app) const
  {
    auto result = preconditionable::template operator()<
      lspg_state_t, fom_t>(romY, app);
    app.applyPreconditioner(*yFom_.data(), *result.data());
    return result;
  }

  template <
      typename lspg_state_t,
      typename lspg_residual_t,
      typename fom_t>
  void operator()(const lspg_state_t  & romY,
                  lspg_residual_t & romR,
                  const fom_t   & app) const
  {
    preconditionable::template operator()<
      lspg_state_t, lspg_residual_t, fom_t>(romY, romR, app);
    app.applyPreconditioner(*yFom_.data(), *romR.data());
  }

};//end class



/* overload when decorating a jacobian policy */
template <typename preconditionable>
class Preconditioned<
  preconditionable,
  ::pressio::mpl::enable_if_t<preconditionable::isResidualPolicy_ == false>
  > : public preconditionable{

public:
  using typename preconditionable::apply_jac_return_t;

protected:
  using preconditionable::yFom_;

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
    app.applyPreconditioner(*yFom_.data(), *JJ.data(), t);
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
    app.applyPreconditioner(*yFom_.data(), *odeJJ.data(), t);
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
    auto JJ = preconditionable::template operator()<
      lspg_state_t, app_t>(romY, app);
    app.applyPreconditioner(*yFom_.data(), *JJ.data());
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
    preconditionable::template operator()<
      lspg_state_t, lspg_jac_t, app_t>(romY, romJJ, app);
    app.applyPreconditioner(*yFom_.data(), *romJJ.data());
  }

};//end class

}}} //end namespace pressio::rom::decorator

#endif
