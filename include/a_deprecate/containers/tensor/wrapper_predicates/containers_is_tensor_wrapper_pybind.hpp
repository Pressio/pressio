/*
//@HEADER
// ************************************************************************
//
// containers_is_tensor_wrapper_pybind.hpp
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

#ifndef CONTAINERS_TENSOR_WRAPPER_PREDICATES_CONTAINERS_IS_TENSOR_WRAPPER_PYBIND_HPP_
#define CONTAINERS_TENSOR_WRAPPER_PREDICATES_CONTAINERS_IS_TENSOR_WRAPPER_PYBIND_HPP_

namespace pressio{ namespace containers{ namespace predicates {

//*******************
// rank-1
//*******************
template <typename T, typename enable = void>
struct is_rank1_tensor_wrapper_pybind : std::false_type {};

template <typename wrapped_t>
struct is_rank1_tensor_wrapper_pybind<
  ::pressio::containers::Tensor<1, wrapped_t>,
  ::pressio::mpl::enable_if_t<
    is_fstyle_array_pybind11<wrapped_t>::value or
    is_cstyle_array_pybind11<wrapped_t>::value
    >
  > : std::true_type{};

template <typename wrapped_t>
struct is_rank1_tensor_wrapper_pybind<
  const ::pressio::containers::Tensor<1, wrapped_t>,
  ::pressio::mpl::enable_if_t<
    is_fstyle_array_pybind11<wrapped_t>::value or
    is_cstyle_array_pybind11<wrapped_t>::value
    >
  > : std::true_type{};


//*******************
// rank-2
//*******************
// fstyle
template <typename T, typename enable = void>
struct is_fstyle_rank2_tensor_wrapper_pybind : std::false_type {};

template <typename wrapped_t>
struct is_fstyle_rank2_tensor_wrapper_pybind<
  ::pressio::containers::Tensor<2, wrapped_t>,
  ::pressio::mpl::enable_if_t<is_fstyle_array_pybind11<wrapped_t>::value>
  > : std::true_type{};

template <typename wrapped_t>
struct is_fstyle_rank2_tensor_wrapper_pybind<
  const ::pressio::containers::Tensor<2, wrapped_t>,
  ::pressio::mpl::enable_if_t<is_fstyle_array_pybind11<wrapped_t>::value>
  > : std::true_type{};

// cstyle
template <typename T, typename enable = void>
struct is_cstyle_rank2_tensor_wrapper_pybind : std::false_type {};

template <typename wrapped_t>
struct is_cstyle_rank2_tensor_wrapper_pybind<
  ::pressio::containers::Tensor<2, wrapped_t>,
  ::pressio::mpl::enable_if_t<is_cstyle_array_pybind11<wrapped_t>::value>
  > : std::true_type{};

template <typename wrapped_t>
struct is_cstyle_rank2_tensor_wrapper_pybind<
  const ::pressio::containers::Tensor<2, wrapped_t>,
  ::pressio::mpl::enable_if_t<is_cstyle_array_pybind11<wrapped_t>::value>
  > : std::true_type{};

//is rank2 tensor
template <typename T, typename enable = void>
struct is_rank2_tensor_wrapper_pybind : std::false_type {};

template <typename T>
struct is_rank2_tensor_wrapper_pybind<
  T,
  ::pressio::mpl::enable_if_t<
    is_fstyle_rank2_tensor_wrapper_pybind<T>::value or
    is_cstyle_rank2_tensor_wrapper_pybind<T>::value
    >
  > : std::true_type{};


//*******************
// rank-3
//*******************
// fstyle
template <typename T, typename enable = void>
struct is_fstyle_rank3_tensor_wrapper_pybind : std::false_type {};

template <typename wrapped_t>
struct is_fstyle_rank3_tensor_wrapper_pybind<
  ::pressio::containers::Tensor<3, wrapped_t>,
  ::pressio::mpl::enable_if_t<is_fstyle_array_pybind11<wrapped_t>::value>
  > : std::true_type{};

template <typename wrapped_t>
struct is_fstyle_rank3_tensor_wrapper_pybind<
  const ::pressio::containers::Tensor<3, wrapped_t>,
  ::pressio::mpl::enable_if_t<is_fstyle_array_pybind11<wrapped_t>::value>
  > : std::true_type{};

// cstyle
template <typename T, typename enable = void>
struct is_cstyle_rank3_tensor_wrapper_pybind : std::false_type {};

template <typename wrapped_t>
struct is_cstyle_rank3_tensor_wrapper_pybind<
  ::pressio::containers::Tensor<3, wrapped_t>,
  ::pressio::mpl::enable_if_t<is_cstyle_array_pybind11<wrapped_t>::value>
  > : std::true_type{};

template <typename wrapped_t>
struct is_cstyle_rank3_tensor_wrapper_pybind<
  const ::pressio::containers::Tensor<3, wrapped_t>,
  ::pressio::mpl::enable_if_t<is_cstyle_array_pybind11<wrapped_t>::value>
  > : std::true_type{};

//is rank3 tensor
template <typename T, typename enable = void>
struct is_rank3_tensor_wrapper_pybind : std::false_type {};

template <typename T>
struct is_rank3_tensor_wrapper_pybind<
  T,
  ::pressio::mpl::enable_if_t<
    is_fstyle_rank3_tensor_wrapper_pybind<T>::value or
    is_cstyle_rank3_tensor_wrapper_pybind<T>::value
    >
  > : std::true_type{};

//--------------------
// fstyle tensor
//--------------------
template <typename T, typename enable = void>
struct is_fstyle_tensor_wrapper_pybind : std::false_type {};

template <typename T>
struct is_fstyle_tensor_wrapper_pybind<
  T,
  ::pressio::mpl::enable_if_t<
    is_rank1_tensor_wrapper_pybind<T>::value or
    is_fstyle_rank2_tensor_wrapper_pybind<T>::value or
    is_fstyle_rank3_tensor_wrapper_pybind<T>::value
    >
  > : std::true_type{};

//--------------------
// cstyle tensor
//--------------------
template <typename T, typename enable = void>
struct is_cstyle_tensor_wrapper_pybind : std::false_type {};

template <typename T>
struct is_cstyle_tensor_wrapper_pybind<
  T,
  ::pressio::mpl::enable_if_t<
    is_rank1_tensor_wrapper_pybind<T>::value or
    is_cstyle_rank2_tensor_wrapper_pybind<T>::value or
    is_cstyle_rank3_tensor_wrapper_pybind<T>::value
    >
  > : std::true_type{};

//--------------------
// is tensor wrapper
//--------------------
template <typename T, typename enable = void>
struct is_tensor_wrapper_pybind : std::false_type {};

template <typename T>
struct is_tensor_wrapper_pybind<
  T,
  ::pressio::mpl::enable_if_t<
    is_fstyle_tensor_wrapper_pybind<T>::value or
    is_cstyle_tensor_wrapper_pybind<T>::value
    >
  > : std::true_type{};

template<typename T>
using is_tensor_wrapper_pybind11 = is_tensor_wrapper_pybind<T>;

}}}//end namespace pressio::containers::predicates
#endif  // CONTAINERS_TENSOR_WRAPPER_PREDICATES_CONTAINERS_IS_TENSOR_WRAPPER_PYBIND_HPP_
