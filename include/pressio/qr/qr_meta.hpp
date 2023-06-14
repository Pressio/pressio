/*
//@HEADER
// ************************************************************************
//
// qr_meta.hpp
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

#ifndef QR_QR_META_HPP_
#define QR_QR_META_HPP_

namespace pressio{ namespace qr{ namespace meta {

#if 0 // not used (`is_dense` and `is_shared_mem` traits are used only here...)
template <typename T, typename Enable = void>
struct is_legitimate_r_type : std::false_type {};

template <typename T>
struct is_legitimate_r_type<T,
	 ::pressio::mpl::enable_if_t<
  	   ::pressio::Traits<T>::rank == 2 and
	   ::pressio::Traits<T>::is_shared_mem and
	   ::pressio::Traits<T>::is_dense
	   >
      > : std::true_type{};
#endif

template <typename T, typename Q_T, typename Enable = void>
struct is_legitimate_vector_type_for_qr_project : std::false_type {};

template <typename T, typename Q_t>
struct is_legitimate_vector_type_for_qr_project<T, Q_t,
	 ::pressio::mpl::enable_if_t<
  	   ::pressio::Traits<T>::rank == 1 and
	   // the vector type should be from same package as Q
	   (false
#ifdef PRESSIO_ENABLE_TPL_EIGEN
	   or (::pressio::is_native_container_eigen<T>::value and ::pressio::is_native_container_eigen<Q_t>::value)
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
	   or (::pressio::is_native_container_kokkos<T>::value and ::pressio::is_native_container_kokkos<Q_t>::value)
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
	   or (::pressio::is_vector_tpetra<T>::value and ::pressio::is_multi_vector_tpetra<Q_t>::value)
	   or (::pressio::is_vector_epetra<T>::value and ::pressio::is_multi_vector_epetra<Q_t>::value)
	   or (::pressio::is_vector_tpetra_block<T>::value and ::pressio::is_multi_vector_tpetra_block<Q_t>::value)
#endif
	   )
	 >
      > : std::true_type{};


// #if defined PRESSIO_ENABLE_TPL_TRILINOS
// template <typename algo_t, typename enable = void>
// struct is_legitimate_algo_for_epetra_mv : std::false_type {};

// template <typename algo_t>
// struct is_legitimate_algo_for_epetra_mv<algo_t,
// 	 ::pressio::mpl::enable_if_t<
// 	   std::is_same<algo_t, ::pressio::qr::Householder>::value
// 	   or std::is_same<algo_t, ::pressio::qr::TSQR>::value
// 	 >
//       > : std::true_type{};
// #endif

// #if defined PRESSIO_ENABLE_TPL_TRILINOS
// template <typename algo_t, typename enable = void>
// struct is_legitimate_algo_for_tpetra_mv : std::false_type {};

// template <typename algo_t>
// struct is_legitimate_algo_for_tpetra_mv<algo_t,
// 	 ::pressio::mpl::enable_if_t<
// 	   std::is_same<algo_t, ::pressio::qr::Householder>::value
// 	   or std::is_same<algo_t, ::pressio::qr::TSQR>::value
// 	   >
//         > : std::true_type{};
// #endif

}}}//end namespace pressio::qr::meta
#endif  // QR_QR_META_HPP_
