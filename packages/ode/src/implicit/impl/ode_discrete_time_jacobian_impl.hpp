/*
//@HEADER
// ************************************************************************
//
// ode_discrete_time_jacobian_impl.hpp
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

#ifndef ODE_IMPLICIT_IMPL_ODE_DISCRETE_TIME_JACOBIAN_IMPL_HPP_
#define ODE_IMPLICIT_IMPL_ODE_DISCRETE_TIME_JACOBIAN_IMPL_HPP_

namespace pressio{ namespace ode{ namespace impl{

#ifdef PRESSIO_ENABLE_TPL_PYBIND11
template <typename jacobian_type, typename scalar_type>
void discrete_time_jacobian(jacobian_type & jac,
			    const scalar_type & dt,
			    ::pressio::ode::implicitmethods::Euler)
{
  /* this is here for compilation purposes only when only pybind11 is enabled
     otherwise it would not compile.
     This should never be called. if it is, it means something is wrong*/
  assert(1==0);
}
#endif

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename jacobian_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  containers::predicates::is_sparse_matrix_wrapper_eigen<jacobian_type>::value
>
discrete_time_jacobian(jacobian_type & jac,
		       const scalar_type & dt,
		       ::pressio::ode::implicitmethods::Euler)
{
  constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_;
  const auto cf	  = ::pressio::ode::constants::bdf1<scalar_type>::c_f_ * dt;
  ::pressio::ops::scale(jac, cf);
  ::pressio::ops::addToDiagonal(jac, cn);
}

template <typename jacobian_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  containers::predicates::is_dense_matrix_wrapper_eigen<jacobian_type>::value
>
discrete_time_jacobian(jacobian_type & jac,
		       const scalar_type & dt,
		       ::pressio::ode::implicitmethods::Euler)
{
  assert(jac.extent(0) == jac.extent(1));
  constexpr auto cn   = ::pressio::ode::constants::bdf1<scalar_type>::c_n_;
  const auto cf	  = ::pressio::ode::constants::bdf1<scalar_type>::c_f_ * dt;
  ::pressio::ops::scale(jac, cf);
  for (auto i=0; i<jac.extent(0); ++i) jac(i,i) += cn;
}

template <typename jacobian_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  containers::predicates::is_sparse_matrix_wrapper_eigen<jacobian_type>::value
  >
discrete_time_jacobian(jacobian_type & jac,
		       const scalar_type & dt,
		       ::pressio::ode::implicitmethods::BDF2)
{
  constexpr auto cn   = ::pressio::ode::constants::bdf2<scalar_type>::c_n_;
  const auto cf	  = ::pressio::ode::constants::bdf2<scalar_type>::c_f_ * dt;
  using namespace ::pressio::ode::constants;
  ::pressio::ops::scale(jac, cf);
  ::pressio::ops::addToDiagonal(jac, cn);
}

template <typename jacobian_type, typename scalar_type>
::pressio::mpl::enable_if_t<
  containers::predicates::is_dense_matrix_wrapper_eigen<jacobian_type>::value
>
discrete_time_jacobian(jacobian_type & jac,
		       const scalar_type & dt,
		       ::pressio::ode::implicitmethods::BDF2)
{
  assert(jac.extent(0) == jac.extent(1));
  constexpr auto cn   = ::pressio::ode::constants::bdf2<scalar_type>::c_n_;
  const auto cf	  = ::pressio::ode::constants::bdf2<scalar_type>::c_f_ * dt;
  ::pressio::ops::scale(jac, cf);
  for (auto i=0; i<jac.extent(0); ++i) jac(i,i) += cn;
}
#endif

}}}//end namespace pressio::ode::impl
#endif  // ODE_IMPLICIT_IMPL_ODE_DISCRETE_TIME_JACOBIAN_IMPL_HPP_
