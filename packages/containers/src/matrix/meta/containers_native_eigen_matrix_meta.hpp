/*
//@HEADER
// ************************************************************************
//
// containers_native_eigen_matrix_meta.hpp
//                     		      Pressio 
// Copyright 2019 National Technology & Engineering Solutions of Sandia,LLC 
//							      (NTESS)
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

#ifndef CONTAINERS_NATIVE_EIGEN_MATRIX_META_HPP_
#define CONTAINERS_NATIVE_EIGEN_MATRIX_META_HPP_

#include "../../meta/containers_meta_basic.hpp"
#include "../../vector/meta/containers_native_eigen_vector_meta.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_dense_static_matrix_eigen : std::false_type {};

/* T is a dense STATIC eigen matrix if
 * T is not an eigen vector
 * rows and cols are not = Eigen:Dynamic
 */
template<typename T>
struct is_dense_static_matrix_eigen<
  T,
  typename
  std::enable_if<
    !is_vector_eigen<T>::value and
    std::is_same<
      T,
      Eigen::Matrix<typename T::Scalar,
		    T::RowsAtCompileTime,
		    T::ColsAtCompileTime
		    >
      >::value and
    T::RowsAtCompileTime != Eigen::Dynamic and
    T::ColsAtCompileTime != Eigen::Dynamic
    >::type
  > : std::true_type{};
//----------------------------------------------------------------------


template <typename T, typename enable = void>
struct is_dense_dynamic_matrix_eigen : std::false_type {};

/* T is a dense DYNAMIC eigen matrix if
 * is not an eigen vector
 * is not a static dense matrix
 * both # rows and cols are dynamic
*/
template<typename T>
struct is_dense_dynamic_matrix_eigen<
  T,
  typename
  std::enable_if<
    !is_vector_eigen<T>::value and
    !is_dense_static_matrix_eigen<T>::value and
    std::is_same<
      T, Eigen::Matrix<typename T::Scalar,
		       Eigen::Dynamic, Eigen::Dynamic
		       >
      >::value
    >::type
  > : std::true_type{};


/* T is a dense DYNAMIC eigen matrix if
 * is not an eigen vector
 * is not a static dense matrix
 * # of rows is static but cols are dynamic
*/
template<typename T>
struct is_dense_dynamic_matrix_eigen<
  T,
  typename
  std::enable_if<
    !is_vector_eigen<T>::value and
    !is_dense_static_matrix_eigen<T>::value and
    std::is_same<
      T, Eigen::Matrix<typename T::Scalar,
		       T::RowsAtCompileTime,
		       Eigen::Dynamic>
      >::value and
    T::RowsAtCompileTime != Eigen::Dynamic
    >::type
  > : std::true_type{};


/* T is a dense DYNAMIC eigen matrix if
 * is not an eigen vector
 * is not a static dense matrix
 * # of rows is dynamic but cols are static
*/
template<typename T>
struct is_dense_dynamic_matrix_eigen<
  T,
  typename
  std::enable_if<
    !is_vector_eigen<T>::value and
    !is_dense_static_matrix_eigen<T>::value and
    std::is_same<
      T, Eigen::Matrix<typename T::Scalar,
		       Eigen::Dynamic,
		       T::ColsAtCompileTime>
      >::value and
    T::ColsAtCompileTime != Eigen::Dynamic
    >::type
  > : std::true_type{};

//----------------------------------------------------------------------


template <typename T, typename enable = void>
struct is_dense_matrix_eigen : std::false_type {};

template<typename T>
struct is_dense_matrix_eigen<
  T, typename
  std::enable_if<
       is_dense_static_matrix_eigen<T>::value or
       is_dense_dynamic_matrix_eigen<T>::value
    >::type
  > : std::true_type{};


//----------------------------------------------------------------------

template <typename T, typename enable = void>
struct is_sparse_matrix_eigen : std::false_type {};

/*
 * T is an eigen sparse matrix if is
 * not an eigen vector
 * is same type as sparse matrix
*/
template<typename T>
struct is_sparse_matrix_eigen<T, typename
  std::enable_if< !is_vector_eigen<T>::value and
		  /*!is_dense_matrix_eigen<T>::value &&*/
		  std::is_same<T,
			       Eigen::SparseMatrix<
			       typename T::Scalar,
				 T::Options,
			       typename T::StorageIndex>
			      >::value
		  >::type
      > : std::true_type{};

//----------------------------------------------------------------------

template <typename T1, typename T2, typename enable = void>
struct sparse_sharedmem_eigen_same_storage : std::false_type{};

template <typename T1, typename T2>
struct sparse_sharedmem_eigen_same_storage<
  T1, T2, typename
  std::enable_if<
	    (T1::is_row_major && T2::is_row_major) ||
	    (T1::is_col_major && T2::is_col_major)
	    >::type
  > : std::true_type{};

//----------------------------------------------------------------------

}}}//end namespace pressio::containers::meta
#endif
