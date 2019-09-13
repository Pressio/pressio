/*
//@HEADER
// ************************************************************************
//
// ode_implicit_jacobian_standard_policy.hpp
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

#ifndef ODE_POLICIES_STANDARD_IMPLICIT_JACOBIAN_STANDARD_POLICY_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_JACOBIAN_STANDARD_POLICY_HPP_

#include "../../../ode_fwd.hpp"
#include "../base/ode_jacobian_policy_base.hpp"
#include "../../ode_jacobian_impl.hpp"

namespace pressio{ namespace ode{ namespace policy{

/*
 * state and jacobian types are containers wrappers
 * both are wrappers from containers
 */
template<
  typename state_type,
  typename system_type,
  typename jacobian_type
  >
class ImplicitJacobianStandardPolicy<
  state_type, system_type, jacobian_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::ode::meta::is_legitimate_implicit_state_type<state_type>::value and
    ::pressio::ode::meta::is_legitimate_jacobian_type<jacobian_type>::value and
    containers::meta::is_wrapper<state_type>::value and
    containers::meta::is_wrapper<jacobian_type>::value
    >
  > : public JacobianPolicyBase<ImplicitJacobianStandardPolicy<
    state_type, system_type, jacobian_type> >
{

  using this_t = ImplicitJacobianStandardPolicy<state_type, system_type, jacobian_type>;
  friend JacobianPolicyBase<this_t>;

public:
  ImplicitJacobianStandardPolicy() = default;
  ~ImplicitJacobianStandardPolicy() = default;

public:
  template <
    ode::ImplicitEnum method, typename scalar_t
  >
  void operator()(const state_type & y,
		  jacobian_type & J,
		  const system_type & model,
		  scalar_t t,
		  scalar_t dt)const
  {
    model.jacobian( *y.data(), t, *J.data());
    ::pressio::ode::impl::time_discrete_jacobian<method>(J, dt);
  }

  template <
    ode::ImplicitEnum method, typename scalar_t
    >
  jacobian_type operator()(const state_type & y,
  			   const system_type & model,
  			   scalar_t t,
  			   scalar_t dt)const
  {
    auto nJJ = model.jacobian(*y.data(), t);
    jacobian_type JJ(nJJ);
    ::pressio::ode::impl::time_discrete_jacobian<method>(JJ, dt);
    return JJ;
  }

};//end class



// #ifdef HAVE_PYBIND11
// /*
//  * state_type and jacobian_type = pybind11::array_t
//  */
// template<typename state_type,
// 	 typename system_type,
// 	 typename jacobian_type>
// class ImplicitJacobianStandardPolicy<
//   state_type, system_type, jacobian_type,
//   ::pressio::mpl::enable_if_t<
//     ::pressio::ode::meta::is_legitimate_implicit_state_type<state_type>::value and
//     ::pressio::ode::meta::is_legitimate_jacobian_type<jacobian_type>::value and
//     mpl::is_same<system_type, pybind11::object >::value and
//     containers::meta::is_array_pybind11<state_type>::value and
//     containers::meta::is_array_pybind11<jacobian_type>::value
//     >
//   > : public JacobianPolicyBase<ImplicitJacobianStandardPolicy<
//     state_type, system_type, jacobian_type> >{

//   using this_t = ImplicitJacobianStandardPolicy<state_type, system_type, jacobian_type>;
//   friend JacobianPolicyBase<this_t>;

// public:
//   ImplicitJacobianStandardPolicy() = default;
//   ~ImplicitJacobianStandardPolicy() = default;

// public:
//   template <
//     ode::ImplicitEnum method, typename scalar_t
//   >
//   void operator()(const state_type & y,
// 		  jacobian_type & J,
// 		  const system_type & model,
// 		  scalar_t t,
// 		  scalar_t dt)const
//   {
//     throw std::runtime_error(" ImplicitJacobianStandardPolicy for pybind11 missing ");
//     // model.jacobian( *y.data(), *J.data(), t);
//     // ::pressio::ode::impl::time_discrete_jacobian<method>(J, dt);
//   }

//   template <
//     ode::ImplicitEnum method, typename scalar_t
//     >
//   jacobian_type operator()(const state_type & y,
//   			   const system_type & model,
//   			   scalar_t t,
//   			   scalar_t dt)const
//   {
//     auto nJJ = model.attr("jacobian")(y, t);
//     throw std::runtime_error(" ImplicitJacobianStandardPolicy for pybind11 missing ");

//     // auto nJJ = model.jacobian(*y.data(), t);
//     // jacobian_type JJ(nJJ);
//     // ::pressio::ode::impl::time_discrete_jacobian<method>(JJ, dt);
//     return nJJ;
//   }
// };//end class
// #endif

}}}//end namespace pressio::ode::policy
#endif
