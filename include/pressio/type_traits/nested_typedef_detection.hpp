/*
//@HEADER
// ************************************************************************
//
// nested_typedef_detection.hpp
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

#ifndef TYPE_TRAITS_NESTED_TYPEDEF_DETECTION_HPP_
#define TYPE_TRAITS_NESTED_TYPEDEF_DETECTION_HPP_

namespace pressio{

template <typename T, typename enable = void>
struct has_state_typedef : std::false_type{};

template <typename T>
struct has_state_typedef<
  T,
  mpl::enable_if_t< !std::is_void<typename T::state_type>::value >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_velocity_typedef : std::false_type{};

template <typename T>
struct has_velocity_typedef<
  T,
  mpl::enable_if_t< !std::is_void<typename T::velocity_type>::value >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_residual_typedef : std::false_type{};

template <typename T>
struct has_residual_typedef<
  T,
  mpl::enable_if_t< !std::is_void<typename T::residual_type>::value >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_scalar_typedef : std::false_type{};

template <typename T>
struct has_scalar_typedef<
  T,
  mpl::enable_if_t< !std::is_void<typename T::scalar_type>::value >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_communicator_typedef : std::false_type{};

template <typename T>
struct has_communicator_typedef<
  T,
  mpl::enable_if_t< !std::is_void<typename T::communicator_type>::value >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_ordinal_typedef : std::false_type{};

template <typename T>
struct has_ordinal_typedef<
  T,
  mpl::enable_if_t< !std::is_void<typename T::ordinal_type>::value >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_matrix_typedef : std::false_type{};

template <typename T>
struct has_matrix_typedef<
  T,
  mpl::enable_if_t< !std::is_void<typename T::matrix_type>::value >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_local_ordinal_typedef : std::false_type{};

template <typename T>
struct has_local_ordinal_typedef<
  T,
  mpl::enable_if_t< !std::is_void<typename T::local_ordinal_type>::value >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_jacobian_typedef : std::false_type{};

template <typename T>
struct has_jacobian_typedef<
  T,
  mpl::enable_if_t< !std::is_void<typename T::jacobian_type>::value >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_hessian_typedef : std::false_type{};

template <typename T>
struct has_hessian_typedef<
  T,
  mpl::enable_if_t< !std::is_void<typename T::hessian_type>::value >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_gradient_typedef : std::false_type{};

template <typename T>
struct has_gradient_typedef<
  T,
  mpl::enable_if_t< !std::is_void<typename T::gradient_type>::value >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_global_ordinal_typedef : std::false_type{};

template <typename T>
struct has_global_ordinal_typedef<
  T,
  mpl::enable_if_t< !std::is_void<typename T::global_ordinal_type>::value >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_data_map_typedef : std::false_type{};

template <typename T>
struct has_data_map_typedef<
  T,
  mpl::enable_if_t< !std::is_void<typename T::data_map_type>::value >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_discrete_time_residual_typedef : std::false_type{};

template <typename T>
struct has_discrete_time_residual_typedef<
  T,
  mpl::enable_if_t< !std::is_void<typename T::discrete_time_residual_type>::value >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_discrete_time_jacobian_typedef : std::false_type{};

template <typename T>
struct has_discrete_time_jacobian_typedef<
  T,
  mpl::enable_if_t< !std::is_void<typename T::discrete_time_jacobian_type>::value >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_dense_matrix_typedef : std::false_type{};

template <typename T>
struct has_dense_matrix_typedef<
  T,
  ::pressio::mpl::enable_if_t<
    !std::is_void<typename T::dense_matrix_type>::value
    >
  > : std::true_type{};

template <typename T, typename enable = void>
struct has_fom_state_typedef : std::false_type{};

template <typename T>
struct has_fom_state_typedef<
  T,
  ::pressio::mpl::enable_if_t<
    !std::is_void<typename T::fom_state_type>::value
    >
  > : std::true_type{};

}//end namespace
#endif  // TYPE_TRAITS_NESTED_TYPEDEF_DETECTION_HPP_
