/*
//@HEADER
// ************************************************************************
//
// containers_native_vector_static_asserts.hpp
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

#ifndef CONTAINERS_NATIVE_VECTOR_STATIC_ASSERTS_HPP_
#define CONTAINERS_NATIVE_VECTOR_STATIC_ASSERTS_HPP_

#include "./meta/containers_native_blaze_vector_meta.hpp"
#include "./meta/containers_native_eigen_vector_meta.hpp"
#include "./meta/containers_native_epetra_vector_meta.hpp"
#include "./meta/containers_native_tpetra_vector_meta.hpp"
#include "./meta/containers_native_teuchos_vector_meta.hpp"
#include "./meta/containers_native_kokkos_vector_meta.hpp"
#include "./meta/containers_native_tpetra_block_vector_meta.hpp"


namespace pressio{ namespace containers{

#define STATIC_ASSERT_IS_VECTOR_EIGEN(TYPE) \
  static_assert( containers::meta::is_vector_eigen<TYPE>::value, \
		 "THIS_IS_NOT_A_VECTOR_FROM_EIGEN")
#define STATIC_ASSERT_IS_NOT_VECTOR_EIGEN(TYPE) \
  static_assert( !containers::meta::is_vector_eigen<TYPE>::value, \
		 "THIS_IS_A_VECTOR_FROM_EIGEN")

#define STATIC_ASSERT_IS_VECTOR_STDLIB(TYPE) \
  static_assert( containers::meta::is_vector_stdlib<TYPE>::value, \
		 "THIS_IS_NOT_A_STDLIB_VECTOR")
#define STATIC_ASSERT_IS_NOT_VECTOR_STDLIB(TYPE) \
  static_assert( !containers::meta::is_vector_stdlib<TYPE>::value, \
		 "THIS_IS_A_STDLIB_VECTOR")

#ifdef HAVE_TRILINOS
#define STATIC_ASSERT_IS_VECTOR_EPETRA(TYPE) \
  static_assert( containers::meta::is_vector_epetra<TYPE>::value, \
		 "THIS_IS_NOT_A_VECTOR_EPETRA")
#define STATIC_ASSERT_IS_NOT_VECTOR_EPETRA(TYPE) \
  static_assert( !containers::meta::is_vector_epetra<TYPE>::value, \
		 "THIS_IS_A_VECTOR_EPETRA")

#define STATIC_ASSERT_IS_VECTOR_TPETRA(TYPE) \
  static_assert( containers::meta::is_vector_tpetra<TYPE>::value, \
		 "THIS_IS_NOT_A_VECTOR_TPETRA")
#define STATIC_ASSERT_IS_NOT_VECTOR_TPETRA(TYPE) \
  static_assert( !containers::meta::is_vector_tpetra<TYPE>::value, \
		 "THIS_IS_A_VECTOR_TPETRA")

#define STATIC_ASSERT_IS_VECTOR_TPETRA_BLOCK(TYPE) \
  static_assert( containers::meta::is_vector_tpetra_block<TYPE>::value, \
		 "THIS_IS_NOT_A_VECTOR_TPETRA_BLOCK")
#define STATIC_ASSERT_IS_NOT_VECTOR_TPETRA_BLOCK(TYPE) \
  static_assert( !containers::meta::is_vector_tpetra_block<TYPE>::value, \
		 "THIS_IS_A_VECTOR_TPETRA_BLOCK")
#endif


#ifdef HAVE_KOKKOS  
#define STATIC_ASSERT_IS_VECTOR_KOKKOS(TYPE) \
  static_assert( containers::meta::is_vector_kokkos<TYPE>::value, \
		 "THIS_IS_NOT_A_VECTOR_KOKKOS")
#define STATIC_ASSERT_IS_NOT_VECTOR_KOKKOS(TYPE)\
  static_assert( !containers::meta::is_vector_kokkos<TYPE>::value, \
		 "THIS_IS_A_VECTOR_KOKKOS")
#endif


#ifdef HAVE_BLAZE
#define STATIC_ASSERT_IS_VECTOR_BLAZE(TYPE)				\
  static_assert( containers::meta::is_static_vector_blaze<TYPE>::value ||	\
                 containers::meta::is_dynamic_vector_blaze<TYPE>::value,	\
		 "THIS_IS_NOT_A_VECTOR_BLAZE")
#define STATIC_ASSERT_IS_NOT_VECTOR_BLAZE(TYPE)				\
  static_assert( !containers::meta::is_static_vector_blaze<TYPE>::value &&	\
                 !containers::meta::is_dynamic_vector_blaze<TYPE>::value,	\
		 "THIS_IS_A_VECTOR_BLAZE")
#endif

#ifdef HAVE_ARMADILLO
#define STATIC_ASSERT_IS_VECTOR_ARMADILLO(TYPE)			\
  static_assert( containers::meta::is_vector_armadillo<TYPE>::value,	\
		 "THIS_IS_NOT_A_VECTOR_ARMADILLO")
#define STATIC_ASSERT_IS_NOT_VECTOR_ARMADILLO(TYPE)		\
  static_assert( !containers::meta::is_vector_armadillo<TYPE>::value,	\
		 "THIS_IS_A_VECTOR_ARMADILLO")
#endif

}}//end namespace pressio::containers
#endif
