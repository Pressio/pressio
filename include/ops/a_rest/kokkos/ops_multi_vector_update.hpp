/*
//@HEADER
// ************************************************************************
//
// ops_multi_vector_update.hpp
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

#ifndef OPS_KOKKOS_OPS_MULTI_VECTOR_UPDATE_HPP_
#define OPS_KOKKOS_OPS_MULTI_VECTOR_UPDATE_HPP_

#include<KokkosBlas1_axpby.hpp>

namespace pressio{ namespace ops{

//----------------------------------------------------------------------
// overloads for computing: MV = a * MV + b * MV1
// where MV is an kokkos multivector wrapper
//----------------------------------------------------------------------
template<typename T1, typename T2, typename scalar_t>
::pressio::mpl::enable_if_t<
  containers::predicates::is_multi_vector_wrapper_kokkos<T1>::value and
  containers::predicates::is_multi_vector_wrapper_kokkos<T2>::value
  >
update(T1 & mv,
	  const scalar_t &a,
	  const T2 & mv1,
	  const scalar_t &b)
{
  /* make sure we don't pass const objects to be modified.
     In kokkos it is legal to modify const views, not for pressio wrappers. */
  static_assert
    (!std::is_const<T1>::value,
     "cannot modify a const-qualified wrapper of a Kokkos view");
  static_assert
    (containers::predicates::are_scalar_compatible<T1, T2>::value,
     "Types are not scalar compatible");
  static_assert
    (::pressio::containers::predicates::have_matching_execution_space<T1, T2>::value,
     "operands need to have same execution space" );

  ::KokkosBlas::axpby(b, *mv1.data(), a, *mv.data());
}

template<typename T1, typename T2, typename scalar_t>
::pressio::mpl::enable_if_t<
  containers::predicates::is_multi_vector_wrapper_kokkos<T1>::value and
  containers::predicates::is_multi_vector_wrapper_kokkos<T2>::value
  >
update(T1 & mv,
	  const T2 & mv1,
	  const scalar_t & b)
{
  /* make sure we don't pass const objects to be modified.
     In kokkos it is legal to modify const views, not for pressio wrappers. */
  static_assert
    (!std::is_const<T1>::value,
     "cannot modify a const-qualified wrapper of a Kokkos view");
  static_assert
    (containers::predicates::are_scalar_compatible<T1, T2>::value,
     "Types are not scalar compatible");
  static_assert
    (::pressio::containers::predicates::have_matching_execution_space<T1, T2>::value,
     "operands need to have same execution space" );

  constexpr auto zero = ::pressio::utils::constants<scalar_t>::zero();
  ::KokkosBlas::axpby(b, *mv1.data(), zero, *mv.data());
}

}}//end namespace pressio::ops
#endif  // OPS_KOKKOS_OPS_MULTI_VECTOR_UPDATE_HPP_
