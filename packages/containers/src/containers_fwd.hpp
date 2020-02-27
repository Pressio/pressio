/*
//@HEADER
// ************************************************************************
//
// containers_fwd.hpp
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

#ifndef CONTAINERS_FORWARD_DECLARATIONS_HPP_
#define CONTAINERS_FORWARD_DECLARATIONS_HPP_

namespace pressio{

struct transpose{};
struct nontranspose{};

struct view{};


namespace containers{

namespace details {
template<typename T, typename enable = void>
struct traits;

template<typename T>
struct traits<const T> : traits<T> {};
}//end namespace containers::details

template<typename derived_type>
class ContainerBase;

template<typename derived_type>
class ContainerDistributedBase;
template<typename derived_type>
class ContainerSharedMemBase;

template<typename derived_type>
class MatrixDistributedBase;
template<typename derived_type>
class MatrixSharedMemBase;

template<typename derived_type>
class MultiVectorDistributedBase;
template<typename derived_type>
class MultiVectorSharedMemBase;

template<typename derived_type>
class VectorDistributedBase;
template<typename derived_type>
class VectorSharedMemBase;


template <
  typename wrapped_type,
  typename Enable = void>
class Vector;

template <
  typename wrapped_type,
  typename Enable = void>
class MultiVector;

template <
  typename wrapped_type,
  typename Enable = void>
class Matrix;


namespace expressions{

template <typename mat_t, typename enable = void>
struct SubspanExpr;

template <typename vec_t, typename enable = void>
struct SpanExpr;
}

}} // end namespace pressio::containers
#endif
