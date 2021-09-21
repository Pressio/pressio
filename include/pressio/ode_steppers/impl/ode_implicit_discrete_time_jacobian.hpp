/*
//@HEADER
// ************************************************************************
//
// ode_implicit_discrete_time_jacobian.hpp
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

/*
  BDF1: J(y_n+1) = I - dt*df_n+1/dy_n+1
  on input jac contains  df_n+1/dy_n+1
  on output, jac contains the time-discrete jacobian
*/
template <typename JacobianType, typename ScalarType>
void discrete_time_jacobian(JacobianType & jac,
                            const ScalarType & dt,
                            ::pressio::ode::BDF1)
{
  constexpr auto cnp1   = ::pressio::ode::constants::bdf1<ScalarType>::c_np1_;
  const auto cf   = ::pressio::ode::constants::bdf1<ScalarType>::c_f_ * dt;
  ::pressio::ops::scale(jac, cf);
  ::pressio::ops::add_to_diagonal(jac, cnp1);
}

/*
  BDF2: J(y_n+1) = I - (2/3)*dt*df_n+1/dy_n+1
  - on input jac contains  df_n+1/dy_n+1
  - on output, jac contains the time-discrete jacobian
*/
template <typename JacobianType, typename ScalarType>
void discrete_time_jacobian(JacobianType & jac,
			    const ScalarType & dt,
			    ::pressio::ode::BDF2)
{
  constexpr auto cnp1   = ::pressio::ode::constants::bdf2<ScalarType>::c_np1_;
  const auto cf   = ::pressio::ode::constants::bdf2<ScalarType>::c_f_ * dt;
  using namespace ::pressio::ode::constants;
  ::pressio::ops::scale(jac, cf);
  ::pressio::ops::add_to_diagonal(jac, cnp1);
}

/*
  CRANK NICOLSON: J(y_n+1) = I - 0.5*dt*df_n+1/dy_n+1
  - on input jac contains  df_n+1/dy_n+1
  - on output, jac contains the time-discrete jacobian
*/
template <typename JacobianType, typename ScalarType>
void discrete_time_jacobian(JacobianType & jac,
			    const ScalarType & dt,
			    ::pressio::ode::CrankNicolson)
{
  using cnst = ::pressio::ode::constants::cranknicolson<ScalarType>;
  constexpr auto cnp1  = cnst::c_np1_;
  const auto cf = cnst::c_fnp1_ * dt;
  ::pressio::ops::scale(jac, cf);
  ::pressio::ops::add_to_diagonal(jac, cnp1);
}

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// /*
//   BDF1 time-discrete jacobian:
//   J(y_n+1) = I - dt*df_n+1/dy_n+1

//   - on input jac contains  df_n+1/dy_n+1
//   - on output, jac contains the time-discrete jacobian
// */
// template <typename T, typename ScalarType>
// ::pressio::mpl::enable_if_t<
//   ::pressio::containers::predicates::is_array_pybind<T>::value
// >
// discrete_time_jacobian(::pressio::containers::Tensor<2, T> & jac,
//           const ScalarType & dt,
//           ::pressio::ode::BDF1)
// {
//   assert(::pressio::ops::extent(jac,0) == ::pressio::ops::extent(jac,1));
//   constexpr auto cnp1   = ::pressio::ode::constants::bdf1<ScalarType>::c_np1_;
//   const auto cf    = ::pressio::ode::constants::bdf1<ScalarType>::c_f_ * dt;
//   ::pressio::ops::scale(jac, cf);
//   for (auto i=0; i<::pressio::ops::extent(jac,0); ++i) jac(i,i) += cnp1;
// }

// /*
//   BDF2 time-discrete jacobian:
//   J(y_n+1) = I - (2/3)*dt*df_n+1/dy_n+1

//   - on input jac contains  df_n+1/dy_n+1
//   - on output, jac contains the time-discrete jacobian
// */
// template <typename T, typename ScalarType>
// ::pressio::mpl::enable_if_t<
//   ::pressio::containers::predicates::is_array_pybind<T>::value
// >
// discrete_time_jacobian(::pressio::containers::Tensor<2, T> & jac,
//           const ScalarType & dt,
//           ::pressio::ode::BDF2)
// {
//   assert(::pressio::ops::extent(jac,0) == ::pressio::ops::extent(jac,1));
//   constexpr auto cnp1   = ::pressio::ode::constants::bdf2<ScalarType>::c_np1_;
//   const auto cf    = ::pressio::ode::constants::bdf2<ScalarType>::c_f_ * dt;
//   ::pressio::ops::scale(jac, cf);
//   for (auto i=0; i<::pressio::ops::extent(jac,0); ++i) jac(i,i) += cnp1;
// }
// #endif

}}}//end namespace pressio::ode::impl
#endif  // ODE_IMPLICIT_IMPL_ODE_DISCRETE_TIME_JACOBIAN_IMPL_HPP_
