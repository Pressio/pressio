/*
//@HEADER
// ************************************************************************
//
// ops_extent.hpp
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

#ifndef OPS_EIGEN_OPS_EXTENT_HPP_
#define OPS_EIGEN_OPS_EXTENT_HPP_

namespace pressio{ namespace ops{

template<class T, class IndexType>
mpl::enable_if_t<
  ::pressio::is_vector_eigen<T>::value,
  decltype(std::declval<const T>().size())
  >
extent(const T & objectIn, const IndexType i)
{
  assert(i==0);
  return objectIn.size();
}

template<class T, class IndexType>
mpl::enable_if_t<
  ::pressio::is_dense_matrix_eigen<T>::value or
  ::pressio::is_sparse_matrix_eigen<T>::value,
  decltype(std::declval<const T>().size())
  >
extent(const T & objectIn, const IndexType i)
{
  if (i==0){
    return objectIn.rows();
  }
  else{
    return objectIn.cols();
  }
}

template<class T, class IndexType>
mpl::enable_if_t<
  ::pressio::is_expression<T>::value and
  ::pressio::Traits<T>::package_identifier == PackageIdentifier::Eigen,
  decltype(std::declval<const T>().extent(0))
  >
extent(const T & objectIn, const IndexType i)
{
  return objectIn.extent(i);
}

}}
#endif  // OPS_EIGEN_OPS_EXTENT_HPP_
