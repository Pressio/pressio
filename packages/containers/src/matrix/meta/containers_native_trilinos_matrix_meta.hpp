/*
//@HEADER
// ************************************************************************
//
// containers_native_trilinos_matrix_meta.hpp
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

#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_NATIVE_TRILINOS_MATRIX_META_HPP_
#define CONTAINERS_NATIVE_TRILINOS_MATRIX_META_HPP_

#include "../../meta/containers_meta_basic.hpp"
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_LocalMap.h>
#include <Tpetra_CrsMatrix_decl.hpp>
#include "Teuchos_SerialDenseMatrix.hpp"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_sparse_matrix_tpetra : std::false_type {};

template <typename T>
struct is_sparse_matrix_tpetra<T,
      typename
      std::enable_if<
  std::is_same<T,
         Tpetra::CrsMatrix<
           typename T::impl_scalar_type,
           typename T::local_ordinal_type,
           typename T::global_ordinal_type,
           typename T::node_type
           >
         >::value
  >::type
      > : std::true_type{};
//--------------------------------------------

template <typename T, typename enable = void>
struct is_sparse_matrix_epetra
  : std::false_type {};

template<typename T>
struct is_sparse_matrix_epetra<T,
    typename std::enable_if<
      std::is_same<T, Epetra_CrsMatrix>::value
      >::type >
  : std::true_type{};

//-------------------------------------------------

template <typename T, typename enable = void>
struct is_dense_matrix_epetra
  : std::false_type {};

template<typename T>
struct is_dense_matrix_epetra<T,
    typename std::enable_if<
      std::is_same<T, Epetra_MultiVector>::value
      >::type >
  : std::true_type{};
//-------------------------------------------------

template <typename T, typename enable = void>
struct is_dense_matrix_teuchos : std::false_type {};

template <typename T>
struct is_dense_matrix_teuchos<T,
    typename std::enable_if<
	std::is_same<T,
	  Teuchos::SerialDenseMatrix<typename T::ordinalType,
				     typename T::scalarType>
	  >::value
	>::type
      > : std::true_type{};
//-------------------------------------------------

template <typename T, typename enable = void>
struct is_dense_matrix_teuchos_rcp : std::false_type {};

template <typename T>
struct is_dense_matrix_teuchos_rcp<T,
    typename std::enable_if<
      is_teuchos_rcp<T>::value and
      is_dense_matrix_teuchos<typename T::element_type>::value
	>::type
      > : std::true_type{};


}}}//end namespace pressio::containers::meta
#endif
#endif
