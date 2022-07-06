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

#ifndef OPS_TPETRA_BLOCK_OPS_POW_HPP_
#define OPS_TPETRA_BLOCK_OPS_POW_HPP_

namespace pressio{ namespace ops{

// y = |x|^exponent, expo>0
template <typename T1, typename T2>
::pressio::mpl::enable_if_t<
  ::pressio::is_vector_tpetra_block<T1>::value and
  ::pressio::is_vector_tpetra_block<T2>::value
  >
abs_pow(T1 & y,
	const T2 & x,
	const typename ::pressio::Traits<T1>::scalar_type & exponent)
{

  using sc_t = typename ::pressio::Traits<T1>::scalar_type;
  using ord_t = typename ::pressio::Traits<T1>::ordinal_type; // local_ordinal_type

  assert(extent(x,0) == extent(y,0));
  assert(exponent > ::pressio::utils::Constants<sc_t>::zero());
  if (exponent < ::pressio::utils::Constants<sc_t>::zero()){
    throw std::runtime_error("this overload is only for exponent > 0");
  }

  auto y_tp = y.getVectorView();
  // I have to constcast here because for block vector getVectorView is non-const
  const auto x_tp = const_cast<T2 &>(x).getVectorView();
  const auto y_kv = y_tp.getLocalViewDevice(Tpetra::Access::OverwriteAllStruct());
  const auto x_kv = x_tp.getLocalViewDevice(Tpetra::Access::ReadOnlyStruct());
  // NOTE that we need the local length of the tpetra view NOT the block
  Kokkos::parallel_for(y_tp.getLocalLength(),
		       KOKKOS_LAMBDA (const ord_t& i){
			 using std::pow;
			 using std::abs;
			 y_kv(i,0) = pow( abs(x_kv(i,0)), exponent);
		       });
}

// y = |x|^exponent, expo<0
template <typename T1, typename T2>
::pressio::mpl::enable_if_t<
  ::pressio::is_vector_tpetra_block<T1>::value and
  ::pressio::is_vector_tpetra_block<T2>::value
  >
abs_pow(T1 & y,
	const T2 & x,
	const typename ::pressio::Traits<T1>::scalar_type & exponent,
	const typename ::pressio::Traits<T1>::scalar_type & eps)
{

  using sc_t = typename ::pressio::Traits<T1>::scalar_type;
  using ord_t = typename ::pressio::Traits<T1>::ordinal_type; // local_ordinal_type

  assert(extent(x,0) == extent(y,0));
  assert(exponent < ::pressio::utils::Constants<sc_t>::zero());
  if (exponent > ::pressio::utils::Constants<sc_t>::zero()){
    throw std::runtime_error("this overload is only for exponent < 0");
  }

  auto y_tp = y.getVectorView();
  // I have to constcast here because for block vector getVectorView is non-const
  const auto x_tp = const_cast<T2 &>(x).getVectorView();
  const auto y_kv = y_tp.getLocalViewDevice(Tpetra::Access::OverwriteAllStruct());
  const auto x_kv = x_tp.getLocalViewDevice(Tpetra::Access::ReadOnlyStruct());

  constexpr auto one = ::pressio::utils::Constants<sc_t>::one();
  const auto expo = -exponent;
  // NOTE that we need the local length of the tpetra view NOT the block
  Kokkos::parallel_for(y_tp.getLocalLength(),
		       KOKKOS_LAMBDA (const ord_t& i){
			 using std::pow;
			 using std::abs;
			 using std::max;
			 y_kv(i,0) = one/max(eps, pow(abs(x_kv(i,0)), expo));
		       });
}

template <typename T>
::pressio::mpl::enable_if_t<
  ::pressio::is_vector_tpetra_block<T>::value
  >
pow(T & x,
    const typename ::pressio::Traits<T>::scalar_type & exponent)
{
  using ord_t = typename ::pressio::Traits<T>::ordinal_type; // local_ordinal_type

  auto x_tpetraview = x.getVectorView();
  auto x_kv = x_tpetraview.getLocalViewDevice(Tpetra::Access::ReadWriteStruct());

  // NOTE that we need the local length of the tpetra view NOT the block
  Kokkos::parallel_for(x_tpetraview.getLocalLength(),
		       KOKKOS_LAMBDA (const ord_t& i){
			 using std::pow;
			 x_kv(i,0) = pow(x_kv(i,0), exponent);
		       });
}

}}//end namespace pressio::ops
#endif  // OPS_TPETRA_BLOCK_OPS_POW_HPP_
