/*
//@HEADER
// ************************************************************************
//
// solvers_gn_hessian_gradient_operators_hg_api.hpp
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

#ifndef SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_HG_API_HPP_
#define SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_HG_API_HPP_

namespace pressio{ namespace solvers{ namespace nonlinear{ namespace impl{

template <typename h_t, typename g_t>
class HessianGradientOperatorsHGApi
{
  using sc_t = typename ::pressio::containers::details::traits<h_t>::scalar_t;
  g_t g_;
  h_t H_;

public:
  HessianGradientOperatorsHGApi() = delete;
  HessianGradientOperatorsHGApi(HessianGradientOperatorsHGApi const &) = default;
  HessianGradientOperatorsHGApi & operator=(HessianGradientOperatorsHGApi const &) = default;
  HessianGradientOperatorsHGApi(HessianGradientOperatorsHGApi &&) = default;
  HessianGradientOperatorsHGApi & operator=(HessianGradientOperatorsHGApi &&) = default;
  ~HessianGradientOperatorsHGApi() = default;

  template <
   typename system_t, typename state_t,
    mpl::enable_if_t<
      pressio::solvers::constraints::system_hessian_gradient<system_t>::value or
      pressio::solvers::constraints::system_fused_hessian_gradient<system_t>::value,
      int
     > = 0
  >
  HessianGradientOperatorsHGApi(const system_t & system, const state_t & state)
    : g_( system.createGradient() ),
      H_( system.createHessian() )
  {
    ::pressio::ops::set_zero(g_);
    ::pressio::ops::set_zero(H_);
  }

public:
  void resetForNewCall()		{ /* no op */ }
  h_t & hessianRef()			{ return H_; }
  g_t & gradientRef()			{ return g_; }
  const h_t & hessianCRef() const	{ return H_; }
  const g_t & gradientCRef() const	{ return g_; }

  sc_t getParameter(std::string key) const {
    throw std::runtime_error("GN HessGrad operators does not have parameters");
    return {};
  }

  template <typename T>
  void setParameter(std::string key, T value) {
    throw std::runtime_error("GN HessGrad operators do not have parameters");
  }

  template< typename system_t, typename state_t>
  void residualNorm(const system_t & system,
		    const state_t & state,
		    sc_t & residualNorm) const
  {
    system.residualNorm(state, ::pressio::Norm::L2, residualNorm);
  }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::constraints::system_hessian_gradient<system_t>::value
    >
  computeOperators(const system_t & sys,
		   const state_t & state,
		   sc_t & residualNorm,
		   bool recomputeSystemJacobian = true)
  {
    if (recomputeSystemJacobian){
      sys.hessian(state, H_);
    }

    sys.gradient(state, g_, ::pressio::Norm::L2,
		 residualNorm,
		 recomputeSystemJacobian);

    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::constants<sc_t>::negOne());
  }

  template<typename system_t, typename state_t>
  mpl::enable_if_t<
    pressio::solvers::constraints::system_fused_hessian_gradient<system_t>::value
    >
  computeOperators(const system_t & sys,
		   const state_t & state,
		   sc_t & residualNorm,
		   bool recomputeSystemJacobian = true)
  {
    sys.hessianAndGradient(state, H_, g_,
			   ::pressio::Norm::L2,
			   residualNorm, recomputeSystemJacobian);

    // scale because of sign convention
    ::pressio::ops::scale(g_, ::pressio::utils::constants<sc_t>::negOne());
  }
};

}}}}
#endif  // SOLVERS_NONLINEAR_IMPL_OPERATORS_SOLVERS_GN_HESSIAN_GRADIENT_OPERATORS_HG_API_HPP_
