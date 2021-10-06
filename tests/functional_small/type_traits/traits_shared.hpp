/*
//@HEADER
// ************************************************************************
//
// type_traits.hpp
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

/*
  Verify values of common container traits
*/

// -------------------------------------------------

/*
    Verifies traits common for all containers
*/
template <
  typename T,
  pressio::PackageIdentifier pack_id,
  int rank,
  bool is_shared_mem,
  bool is_dynamic,
  typename Scalar,
  typename Ordinal,
  typename SizeType = Ordinal,
  typename ScalarRef = typename std::add_lvalue_reference<
    Scalar
  >::type,
  typename traits = pressio::Traits<T>
>
void test_container_traits()
{
  // ContainersSharedTraits
  static_assert(traits::package_identifier == pack_id);
  static_assert(traits::is_shared_mem == is_shared_mem);
  static_assert(traits::is_distributed == !is_shared_mem);
  static_assert(traits::rank == rank);
  // AllocTrait
  static_assert(traits::is_static == !is_dynamic);
  static_assert(traits::is_dynamic == is_dynamic);
  // ScalarTrait
  testing::StaticAssertTypeEq<typename traits::scalar_type, Scalar>();
  testing::StaticAssertTypeEq<typename traits::reference_type, Scalar &>();
  testing::StaticAssertTypeEq<typename traits::const_reference_type, Scalar const &>();
  // OrdinalTrait
  testing::StaticAssertTypeEq<typename traits::ordinal_type, Ordinal>();
  testing::StaticAssertTypeEq<typename traits::size_type, SizeType>();
}

// -------------------------------------------------

/*
    Verifies traits common for all matrices
*/
template <
  typename T,
  pressio::MatrixIdentifier mtx_id,
  bool is_row_major = true,
  bool is_sparse = false,
  typename traits = pressio::Traits<T>
>
void test_matrix_traits()
{
  static_assert(traits::matrix_identifier == mtx_id);
  static_assert(traits::is_sparse == is_sparse);
  static_assert(traits::is_dense == !is_sparse);
  static_assert(traits::is_row_major == is_row_major);
  static_assert(traits::is_col_major == !is_row_major);
}

// -------------------------------------------------
