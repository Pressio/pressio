/*
//@HEADER
// ************************************************************************
//
// containers_native_multi_vector_static_asserts.hpp
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

#ifndef CONTAINERS_NATIVE_MULTI_VECTOR_STATIC_ASSERTS_HPP_
#define CONTAINERS_NATIVE_MULTI_VECTOR_STATIC_ASSERTS_HPP_

#include "./meta/containers_native_eigen_multi_vector_meta.hpp"
#include "./meta/containers_native_epetra_multi_vector_meta.hpp"
#include "./meta/containers_native_tpetra_multi_vector_meta.hpp"
#include "./meta/containers_native_tpetra_block_multi_vector_meta.hpp"

namespace pressio{ namespace containers{

#ifdef HAVE_TRILINOS
#define STATIC_ASSERT_IS_MULTIVECTOR_EPETRA(TYPE) \
  static_assert( containers::meta::is_multi_vector_epetra<TYPE>::value, \
		 "THIS_IS_NOT_A_MULTIVECTOR_EPETRA")
#define STATIC_ASSERT_IS_NOT_MULTIVECTOR_EPETRA(TYPE) \
  static_assert( !containers::meta::is_multi_vector_epetra<TYPE>::value, \
		 "THIS_IS_A_MULTIVECTOR_EPETRA")

#define STATIC_ASSERT_IS_MULTIVECTOR_TPETRA(TYPE)		  \
  static_assert( containers::meta::is_multi_vector_tpetra<TYPE>::value, \
		 "THIS_IS_NOT_A_MULTIVECTOR_TPETRA")
#define STATIC_ASSERT_IS_NOT_MULTIVECTOR_TPETRA(TYPE) \
  static_assert( !containers::meta::is_multi_vector_tpetra<TYPE>::value, \
		 "THIS_IS_A_MULTIVECTOR_TPETRA")

#define STATIC_ASSERT_IS_MULTIVECTOR_TPETRA_BLOCK(TYPE)		  \
  static_assert( containers::meta::is_multi_vector_tpetra_block<TYPE>::value, \
		 "THIS_IS_NOT_A_MULTIVECTOR_TPETRA_BLOCK")
#define STATIC_ASSERT_IS_NOT_MULTIVECTOR_TPETRA_BLOCK(TYPE) \
  static_assert( !containers::meta::is_multi_vector_tpetra_block<TYPE>::value, \
		 "THIS_IS_A_MULTIVECTOR_TPETRA_BLOCK")

#endif

}}//end namespace pressio::containers
#endif
