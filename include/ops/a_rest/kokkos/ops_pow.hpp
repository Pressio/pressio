/*
//@HEADER
// ************************************************************************
//
// ops_pow.hpp
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

#ifndef OPS_KOKKOS_OPS_POW_HPP_
#define OPS_KOKKOS_OPS_POW_HPP_

namespace pressio{ namespace ops{

template <typename T>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T>::value
  >
pow(T & x,
    const typename ::pressio::containers::details::traits<T>::scalar_t & exponent)
{
  using ord_t = typename ::pressio::containers::details::traits<T>::size_t;
  auto x_kv = *x.data();
  Kokkos::parallel_for(x.extent(0),
		       KOKKOS_LAMBDA (const ord_t& i){
			 using std::pow;
			 x_kv(i) = pow(x_kv(i), exponent);
		       });
}

// y = |x|^exponent, expo>0
template <typename T1, typename T2>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T1>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T2>::value
  >
abs_pow(T1 & y,
	const T2 & x,
	const typename ::pressio::containers::details::traits<T1>::scalar_t & exponent)
{
  static_assert
    (::pressio::containers::predicates::are_scalar_compatible<T1,T2>::value,
     "not scalar compatible");
  using sc_t = typename ::pressio::containers::details::traits<T1>::scalar_t;
  using ord_t = typename ::pressio::containers::details::traits<T1>::size_t;

  assert(x.extent(0) == y.extent(0));
  assert(exponent > ::pressio::utils::constants<sc_t>::zero());
  if (exponent < ::pressio::utils::constants<sc_t>::zero())
    throw std::runtime_error("This overload only supports exponent > 0");

  auto x_kv = *x.data();
  auto y_kv = *y.data();
  Kokkos::parallel_for(x.extent(0),
		       KOKKOS_LAMBDA (const ord_t& i){
			 using std::pow;
			 using std::abs;
			 y_kv(i) = pow( abs(x_kv(i)), exponent);
		       });
}

// y = |x|^exponent, expo<0
template <typename T1, typename T2>
::pressio::mpl::enable_if_t<
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T1>::value and
  ::pressio::ops::constraints::rank1_container_kokkos_with_native_data_access<T2>::value
  >
abs_pow(T1 & y,
	const T2 & x,
	const typename ::pressio::containers::details::traits<T1>::scalar_t & exponent,
	const typename ::pressio::containers::details::traits<T1>::scalar_t & eps)
{
  static_assert
    (::pressio::containers::predicates::are_scalar_compatible<T1,T2>::value,
     "not scalar compatible");
  using sc_t = typename ::pressio::containers::details::traits<T1>::scalar_t;
  using ord_t = typename ::pressio::containers::details::traits<T1>::size_t;

  assert(x.extent(0) == y.extent(0));
  assert(exponent < ::pressio::utils::constants<sc_t>::zero());
  if (exponent > ::pressio::utils::constants<sc_t>::zero())
    throw std::runtime_error("This overload only supports exponent < 0");

  constexpr auto one = ::pressio::utils::constants<sc_t>::one();
  const auto expo = -exponent;
  auto x_kv = *x.data();
  auto y_kv = *y.data();
  Kokkos::parallel_for(x.extent(0),
		       KOKKOS_LAMBDA (const ord_t& i){
			 using std::pow;
			 using std::abs;
			 using std::max;
			 y_kv(i) = one/max(eps, pow(abs(x_kv(i)), expo));
		       });
}

}}//end namespace pressio::ops
#endif  // OPS_KOKKOS_OPS_POW_HPP_
