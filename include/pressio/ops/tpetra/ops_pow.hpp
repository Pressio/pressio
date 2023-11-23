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

#ifndef OPS_TPETRA_OPS_POW_HPP_
#define OPS_TPETRA_OPS_POW_HPP_

namespace pressio{ namespace ops{

// y = |x|^exponent, expo>0
template <typename T1, typename T2>
::pressio::mpl::enable_if_t<
  (::pressio::is_vector_tpetra<T1>::value
   || ::pressio::is_expression_column_acting_on_tpetra<T1>::value)
  && (::pressio::is_vector_tpetra<T2>::value
   || ::pressio::is_expression_column_acting_on_tpetra<T2>::value)
  >
abs_pow(T1 & y_in,
	const T2 & x_in,
	const typename ::pressio::Traits<T1>::scalar_type & exponent)
{

  using sc_t = typename ::pressio::Traits<T1>::scalar_type;

  auto x = impl::get_native(x_in);
  auto y = impl::get_native(y_in);

  assert(x.getGlobalLength() == y.getGlobalLength());
  assert(x.getLocalLength() == y.getLocalLength());
  assert(exponent > ::pressio::utils::Constants<sc_t>::zero());
  if (exponent < ::pressio::utils::Constants<sc_t>::zero()){
    throw std::runtime_error("this overload is only for exponent > 0");
  }

  auto x_kv = x.getLocalViewDevice(Tpetra::Access::ReadOnlyStruct());
  auto y_kv = y.getLocalViewDevice(Tpetra::Access::OverwriteAllStruct());
  Kokkos::parallel_for(x.getLocalLength(),
		       KOKKOS_LAMBDA (const std::size_t& i){
			 using std::pow;
			 using std::abs;
			 y_kv(i,0) = pow( abs(x_kv(i,0)), exponent);
		       });
}

// y = |x|^exponent, expo<0
template <typename T1, typename T2>
::pressio::mpl::enable_if_t<
  (::pressio::is_vector_tpetra<T1>::value
   || ::pressio::is_expression_column_acting_on_tpetra<T1>::value)
  && (::pressio::is_vector_tpetra<T2>::value
   || ::pressio::is_expression_column_acting_on_tpetra<T2>::value)
  >
abs_pow(T1 & y_in,
	const T2 & x_in,
	const typename ::pressio::Traits<T1>::scalar_type & exponent,
	const typename ::pressio::Traits<T1>::scalar_type & eps)
{

  using sc_t = typename ::pressio::Traits<T1>::scalar_type;

  auto x = impl::get_native(x_in);
  auto y = impl::get_native(y_in);

  assert(x.getGlobalLength() == y.getGlobalLength());
  assert(x.getLocalLength() == y.getLocalLength());
  assert(exponent < ::pressio::utils::Constants<sc_t>::zero());
  if (exponent > ::pressio::utils::Constants<sc_t>::zero()){
    throw std::runtime_error("this overload is only for exponent < 0");
  }

  constexpr auto one = ::pressio::utils::Constants<sc_t>::one();
  const auto expo = -exponent;
  auto x_kv = x.getLocalViewDevice(Tpetra::Access::ReadOnlyStruct());
  auto y_kv = y.getLocalViewDevice(Tpetra::Access::OverwriteAllStruct());
  Kokkos::parallel_for(x.getLocalLength(),
		       KOKKOS_LAMBDA (const std::size_t& i){
			 using std::pow;
			 using std::abs;
			 using std::max;
			 y_kv(i,0) = one/max(eps, pow(abs(x_kv(i,0)), expo));
		       });
}

template <typename T>
::pressio::mpl::enable_if_t<
  ::pressio::is_vector_tpetra<T>::value
  || ::pressio::is_expression_column_acting_on_tpetra<T>::value
  >
pow(T & x_in,
    const typename ::pressio::Traits<T>::scalar_type & exponent)
{
  auto x = impl::get_native(x_in);

  auto x_kv = x.getLocalViewDevice(Tpetra::Access::ReadWriteStruct());
  Kokkos::parallel_for(x.getLocalLength(),
		       KOKKOS_LAMBDA (const std::size_t& i){
			 x_kv(i,0) = std::pow(x_kv(i,0), exponent);
		       });
}

}}//end namespace pressio::ops
#endif  // OPS_TPETRA_OPS_POW_HPP_
