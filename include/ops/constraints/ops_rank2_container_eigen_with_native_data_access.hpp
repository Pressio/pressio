/*
//@HEADER
// ************************************************************************
//
// ops_rank2_container_eigen_with_native_data_access.hpp
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

#ifndef OPS_CONSTRAINTS_OPS_RANK2_CONTAINER_EIGEN_WITH_NATIVE_DATA_ACCESS_HPP_
#define OPS_CONSTRAINTS_OPS_RANK2_CONTAINER_EIGEN_WITH_NATIVE_DATA_ACCESS_HPP_

namespace pressio{ namespace ops{ namespace constraints{

template <typename T, typename = void>
struct rank2_container_eigen_with_native_data_access : std::false_type{};

template <typename T>
struct rank2_container_eigen_with_native_data_access<
  T,
  mpl::enable_if_t<
    ::pressio::containers::predicates::is_dense_matrix_wrapper_eigen<T>::value or
    ::pressio::containers::predicates::is_multi_vector_wrapper_eigen<T>::value
    >
  > : std::true_type{};

template <typename T>
struct rank2_container_eigen_with_native_data_access<
  T,
  mpl::enable_if_t<
     ::pressio::containers::predicates::subspan_expression<T>::value and
    T::traits::wrapped_package_identifier ==
    ::pressio::containers::details::WrappedPackageIdentifier::Eigen
    >
  > : std::true_type{};

}}}
#endif  // OPS_CONSTRAINTS_OPS_RANK2_CONTAINER_EIGEN_WITH_NATIVE_DATA_ACCESS_HPP_
