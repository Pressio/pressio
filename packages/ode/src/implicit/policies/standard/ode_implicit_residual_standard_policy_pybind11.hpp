/*
//@HEADER
// ************************************************************************
//
// ode_implicit_residual_standard_policy_pybind11.hpp
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

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#ifndef ODE_POLICIES_STANDARD_IMPLICIT_RESIDUAL_STANDARD_POLICY_PYBIND11_HPP_
#define ODE_POLICIES_STANDARD_IMPLICIT_RESIDUAL_STANDARD_POLICY_PYBIND11_HPP_

#include "../../../ode_fwd.hpp"
#include "../base/ode_implicit_residual_policy_base.hpp"
#include "../../ode_residual_impl.hpp"

namespace pressio{ namespace ode{ namespace policy{

template<
  typename state_type,
  typename system_type,
  typename residual_type
  >
class ImplicitResidualStandardPolicyPybind11<
  state_type, system_type, residual_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::ode::meta::is_legitimate_implicit_state_type<state_type>::value and
    ::pressio::ode::meta::is_legitimate_implicit_residual_type<residual_type>::value and
    mpl::is_same<system_type, pybind11::object >::value and
    containers::meta::is_array_pybind11<state_type>::value and
    containers::meta::is_array_pybind11<residual_type>::value
    >
  >
  : public ImplicitResidualPolicyBase<
  ImplicitResidualStandardPolicyPybind11<state_type, system_type, residual_type>>
{

  using this_t = ImplicitResidualStandardPolicyPybind11<state_type, system_type, residual_type>;
  friend ImplicitResidualPolicyBase<this_t>;

public:
  ImplicitResidualStandardPolicyPybind11() = default;
  ~ImplicitResidualStandardPolicyPybind11() = default;

public:
  template <
    ode::ImplicitEnum method, std::size_t n, typename scalar_type
  >
  void operator()(const state_type & y,
		  residual_type & R,
		  const ::pressio::ode::StatesContainer<state_type, n> & oldYs,
		  const system_type & model,
		  const scalar_type & t,
		  const scalar_type & dt,
		  const types::step_t & step) const
  {
    throw std::runtime_error("ImplicitResidualStandardPolicyPybind11 missing");
    // printf("C++ R address: %p\n", R.data());
    // model.attr("residual2")(y, R, t);
    // ::pressio::ode::impl::time_discrete_residual<method, n>(y, R, oldYs, dt);
  }

  template <
    ode::ImplicitEnum method, std::size_t n, typename scalar_type
    >
  residual_type operator()(const state_type & y,
  			   const ::pressio::ode::StatesContainer<state_type, n> & oldYs,
  			   const system_type & model,
  			   const scalar_type & t,
  			   const scalar_type & dt,
			   const types::step_t & step) const
  {
    throw std::runtime_error("ImplicitResidualStandardPolicyPybind11 missing");
    residual_type nR;// = model.attr("residual1")(y, t);
    // ::pressio::ode::impl::time_discrete_residual<method, n>(y, nR, oldYs, dt);
    return nR;
  }
};//end class

}}}//end namespace pressio::ode::policy
#endif
#endif
