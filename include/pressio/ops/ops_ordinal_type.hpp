/*
//@HEADER
// ************************************************************************
//
// ops_get_native.hpp
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

#ifndef OPS_OPS_ORDINAL_TYPE_HPP_
#define OPS_OPS_ORDINAL_TYPE_HPP_

namespace pressio{ namespace ops{ namespace impl{

template<typename T, class Enable = void> struct OrdinalType       { using type = void; };
template<typename T, class Enable = void> struct LocalOrdinalType  { using type = void; };
template<typename T, class Enable = void> struct GlobalOrdinalType { using type = void; };

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <class T>
struct OrdinalType<
  T,
  std::enable_if_t<
    ::pressio::is_vector_eigen<T>::value
    || ::pressio::is_dense_matrix_eigen<T>::value
    >
  >
{
  using type = typename T::StorageIndex;
};

// following definitions consider native Eigen expression types
template <class VectorType, int Size>
struct OrdinalType<
  Eigen::VectorBlock<VectorType, Size>
  >
: public OrdinalType<VectorType>
{};

template <class MatrixType, int DiagIndex>
struct OrdinalType<
  Eigen::Diagonal<MatrixType, DiagIndex>
  >
: public OrdinalType<MatrixType>
{};

template<class XprType, int BlockRows, int BlockCols, bool InnerPanel>
struct OrdinalType<
  Eigen::Block<XprType, BlockRows, BlockCols, InnerPanel>
  >
: public OrdinalType<XprType>
{};
#endif

#ifdef PRESSIO_ENABLE_TPL_KOKKOS
template <class T>
struct OrdinalType<
  T,
  std::enable_if_t<
    ::pressio::is_vector_kokkos<T>::value
    || ::pressio::is_dense_matrix_kokkos<T>::value
    >
  >
{
  using type = typename T::traits::size_type;
};
#endif

#ifdef PRESSIO_ENABLE_TPL_TRILINOS
template <class T>
struct LocalOrdinalType<
  T,
  std::enable_if_t<
    ::pressio::is_vector_tpetra<T>::value
    || ::pressio::is_vector_tpetra_block<T>::value
    || ::pressio::is_multi_vector_tpetra_block<T>::value
    || ::pressio::is_multi_vector_tpetra<T>::value
    >
 >
{
  using type = typename T::local_ordinal_type;
};

template <class T>
struct GlobalOrdinalType<
  T,
  std::enable_if_t<
    ::pressio::is_vector_tpetra<T>::value
    || ::pressio::is_vector_tpetra_block<T>::value
    || ::pressio::is_multi_vector_tpetra_block<T>::value
    || ::pressio::is_multi_vector_tpetra<T>::value
    >
 >
{
  using type = typename T::global_ordinal_type;
};

template <class T>
struct LocalOrdinalType<
  T,
  std::enable_if_t<
    ::pressio::is_vector_epetra<T>::value
    || ::pressio::is_multi_vector_epetra<T>::value
    >
 >
{
  using type = int;
};

template <class T>
struct GlobalOrdinalType<
  T,
  std::enable_if_t<
    ::pressio::is_vector_epetra<T>::value
    || ::pressio::is_multi_vector_epetra<T>::value
    >
 >
{
  using type = int;
};
#endif

template <class T>
struct OrdinalType<
  T,
  std::enable_if_t<
    ::pressio::is_expression<T>::value
    >
  >
{
private:
  using native_type = typename ::pressio::Traits<T>::native_expr_type;
public:
  using type = typename OrdinalType<native_type>::type;
};

template <class T>
struct LocalOrdinalType<
  T,
  std::enable_if_t<
    ::pressio::is_expression<T>::value
    >
  >
{
private:
  using native_type = typename ::pressio::Traits<T>::native_expr_type;
public:
  using type = typename OrdinalType<native_type>::type;
};

template <class T>
struct GlobalOrdinalType<
  T,
  std::enable_if_t<
    ::pressio::is_expression<T>::value
    >
  >
{
private:
  using native_type = typename ::pressio::Traits<T>::native_expr_type;
public:
  using type = typename OrdinalType<native_type>::type;
};

template <class T>
using ordinal_t = typename OrdinalType<T>::type;
template <class T>
using local_ordinal_t = typename LocalOrdinalType<T>::type;
template <class T>
using global_ordinal_t = typename GlobalOrdinalType<T>::type;

}}}
#endif  // OPS_OPS_ORDINAL_TYPE_HPP_
