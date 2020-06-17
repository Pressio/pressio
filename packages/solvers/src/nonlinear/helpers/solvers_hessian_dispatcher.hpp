/*
//@HEADER
// ************************************************************************
//
// solvers_hessian_dispatcher.hpp
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

#ifndef SOLVERS_IMPL_HESSIAN_DISPATCHER_HPP_
#define SOLVERS_IMPL_HESSIAN_DISPATCHER_HPP_

namespace pressio{ namespace solvers{ namespace iterative{ namespace impl{

template <typename ud_ops_t, typename enable = void>
struct HessianDispatcher;

template <>
struct HessianDispatcher<void>
{

  template <typename J_t, typename result_t>
  void evaluate(const J_t & J, result_t & result) const
  {
    using scalar_t = typename ::pressio::containers::details::traits<J_t>::scalar_t;
    constexpr auto beta  = ::pressio::utils::constants<scalar_t>::zero();
    constexpr auto alpha = ::pressio::utils::constants<scalar_t>::one();
    // here to compute hessian we use the overload taking only J, since we know H = J^T J
    // and that overload is more efficient for doing J^T J since J is the same
    ::pressio::ops::product(::pressio::transpose(), ::pressio::nontranspose(), alpha, J, beta, result);
  }

  template <typename J_t, typename result_t>
  result_t evaluate(const J_t & J) const
  {
    using scalar_t = typename ::pressio::containers::details::traits<J_t>::scalar_t;
    constexpr auto alpha = ::pressio::utils::constants<scalar_t>::one();
    return ::pressio::ops::product<result_t>(::pressio::transpose(),
					     ::pressio::nontranspose(), alpha, J);
  }
};


/*********************
 * user-defined ops
 *********************/
template <typename ud_ops_t>
struct HessianDispatcher<
  ud_ops_t,
  mpl::enable_if_t< !std::is_void<ud_ops_t>::value>
  >
{
private:
  const ud_ops_t * udOps_ = nullptr;

public:
  HessianDispatcher() = delete;
  HessianDispatcher(const ud_ops_t * udOps) : udOps_{udOps}{}

public:
  template < typename J_t, typename result_t, typename _ud_ops_t = ud_ops_t >
  mpl::enable_if_t< !std::is_void<ud_ops_t>::value, result_t>
  evaluate(const J_t & J) const
  {
    using scalar_t	 = typename ::pressio::containers::details::traits<J_t>::scalar_t;
    constexpr auto alpha = ::pressio::utils::constants<scalar_t>::one();
    return udOps_->template product<result_t>(::pressio::transpose(), ::pressio::nontranspose(),
					      alpha, *J.data(), *J.data());
  }

  template <typename J_t, typename result_t, typename _ud_ops_t = ud_ops_t>
  mpl::enable_if_t< !std::is_void<ud_ops_t>::value >
  evaluate(const J_t & J, result_t & result) const
  {
    using scalar_t	 = typename ::pressio::containers::details::traits<J_t>::scalar_t;
    constexpr auto alpha = ::pressio::utils::constants<scalar_t>::one();
    constexpr auto beta  = ::pressio::utils::constants<scalar_t>::zero();
    udOps_->product(::pressio::transpose(), ::pressio::nontranspose(),
			   alpha, *J.data(), *J.data(), beta, result);
  }
};

}}}} //end namespace pressio::solvers::iterative::impl
#endif
