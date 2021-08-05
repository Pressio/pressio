/*
//@HEADER
// ************************************************************************
//
// traits_shared.hpp
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

#ifndef TYPE_TRAITS_TRAITS_SHARED_HPP_
#define TYPE_TRAITS_TRAITS_SHARED_HPP_

namespace pressio{
namespace impl{

template<
  PackageIdentifier PackId,
  bool _is_shared_mem,
  int _rank
  >
struct ContainersSharedTraits
{
  static constexpr PackageIdentifier package_identifier = PackId;
  static constexpr bool is_shared_mem	= _is_shared_mem;
  static constexpr bool is_distributed	= !is_shared_mem;
  static constexpr int rank = _rank;
};

template<
  typename Ordinal,
  typename Size = Ordinal
>
struct OrdinalTrait
{
  using ordinal_type = Ordinal;
  using size_type    = Size;
};

template<
  typename Scalar,
  typename ScalarRef = Scalar&,
  typename ScalarConstRef = typename std::add_const<ScalarRef>::type
>
struct ScalarTrait
{
  using scalar_type    = Scalar;
  using reference_type = ScalarRef;
  using const_reference_type = ScalarConstRef;
};

template<bool _is_static>
struct AllocTrait
{
  static constexpr bool is_static = _is_static;
  static constexpr bool is_dynamic  = !is_static;
};

using StaticAllocTrait = AllocTrait<true>;
using DynamicAllocTrait = AllocTrait<false>;

/// common traits of matrices
template<bool is_sparse_b>
struct MatrixDensityTrait
{
  static constexpr bool is_sparse = is_sparse_b;
  static constexpr bool is_dense  = !is_sparse_b;
};

using DenseMatrixTrait = MatrixDensityTrait<false>;
using SparseMatrixTrait = MatrixDensityTrait<true>;

}
}
#endif  // TYPE_TRAITS_TRAITS_SHARED_HPP_
