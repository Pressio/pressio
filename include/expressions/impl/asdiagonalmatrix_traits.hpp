/*
//@HEADER
// ************************************************************************
//
// containers_asdiagonalmatrix_traits.hpp
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

#ifndef CONTAINERS_EXPRESSIONS_ASDIAGONALMATRIX_CONTAINERS_ASDIAGONALMATRIX_TRAITS_HPP_
#define CONTAINERS_EXPRESSIONS_ASDIAGONALMATRIX_CONTAINERS_ASDIAGONALMATRIX_TRAITS_HPP_

namespace pressio{ namespace expressions{  namespace impl{

#ifdef PRESSIO_ENABLE_TPL_EIGEN
template <typename T>
struct asdiagmatrix_traits<
  AsDiagonalMatrixExpr<T>,
  ::pressio::mpl::enable_if_t<
    ::pressio::is_dynamic_vector_eigen<
      typename std::remove_cv<T>::type
    >::value
    >
  >
  : public containers_shared_traits<PackageIdentifier::Eigen, true, 2>,
    public matrix_shared_traits<false>
{
  static constexpr bool is_static = true;
  static constexpr bool is_dynamic = false;
  using scalar_type  = typename traits<T>::scalar_type;
  using ordinal_type = typename traits<T>::ordinal_type;
  using size_type    = ordinal_type;

  // conditiona ref type because native expression returns by value when object is const
  using reference_type = typename std::conditional<
    std::is_const<T>::value, scalar_type, scalar_type &
  >::type;

  using const_reference_type = typename std::conditional<
    std::is_const<T>::value, scalar_type, scalar_type const &
    >::type;
};
#endif

// #ifdef PRESSIO_ENABLE_TPL_PYBIND11
// template <typename T>
// struct traits<
//   ::pressio::expressions::AsDiagonalMatrixExpr<T>,
//   ::pressio::mpl::enable_if_t<
//     ::pressio::is_rank1_tensor_pybind<T>::value
//     >
//   >
//   : public containers_shared_traits<PackageIdentifier::Pybind, true, 2>,
//     public matrix_shared_traits<false>
// {
//   static constexpr bool is_static = true;
//   static constexpr bool is_dynamic  = !is_static;
//   using scalar_t  = typename traits<T>::scalar_t;
//   using ordinal_t = typename traits<T>::ordinal_t;
//   using size_t    = ordinal_t;

//   // conditional ref type because native expression returns by value when object is const
//   using reference_t = typename traits<T>::reference_t;
//   using const_reference_t = typename traits<T>::const_reference_t;
// };
// #endif

// #ifdef PRESSIO_ENABLE_TPL_TRILINOS
// template <typename T>
// struct traits<
//   ::pressio::expressions::AsDiagonalMatrixExpr<T>,
//   ::pressio::mpl::enable_if_t<
//     ::pressio::is_vector_tpetra<
//     typename std::remove_cv<T>::type
//     >::value
//     >
//   >
//   : public containers_shared_traits<PackageIdentifier::Trilinos, true, 2>,
//     public matrix_shared_traits<false>
// {
//   static constexpr bool is_static = true;
//   static constexpr bool is_dynamic = false;
//   using scalar_t  = typename traits<T>::scalar_t;
//   using local_ordinal_t  = typename traits<T>::local_ordinal_t;
//   using global_ordinal_t = typename traits<T>::global_ordinal_t;
//   using size_t = typename traits<T>::size_t;
// };

// template <typename T>
// struct traits<
//   ::pressio::expressions::AsDiagonalMatrixExpr<T>,
//   ::pressio::mpl::enable_if_t<
//     ::pressio::is_vector_tpetra_block<
//     typename std::remove_cv<T>::type
//     >::value
//     >
//   >
//   : public containers_shared_traits<PackageIdentifier::Trilinos, true, 2>,
//     public matrix_shared_traits<false>
// {
//   static constexpr bool is_static = true;
//   static constexpr bool is_dynamic = false;
//   using wrapped_t = typename traits<T>::wrapped_t;
//   using scalar_t  = typename traits<T>::scalar_t;
//   using local_ordinal_t  = typename traits<T>::local_ordinal_t;
//   using global_ordinal_t = typename traits<T>::global_ordinal_t;
//   using size_t = typename traits<T>::size_t;
// };

// template <typename T>
// struct traits<
//   ::pressio::expressions::AsDiagonalMatrixExpr<T>,
//   ::pressio::mpl::enable_if_t<
//     ::pressio::is_vector_epetra<
//     typename std::remove_cv<T>::type
//     >::value
//     >
//   >
//   : public containers_shared_traits<PackageIdentifier::Trilinos, true, 2>,
//     public matrix_shared_traits<false>
// {

//   static constexpr bool is_static = true;
//   static constexpr bool is_dynamic = false;
//   using scalar_t  = typename traits<T>::scalar_t;
//   using local_ordinal_t  = typename traits<T>::local_ordinal_t;
//   using global_ordinal_t = typename traits<T>::global_ordinal_t;
//   using size_t = typename traits<T>::size_t;
// };
// #endif

}}}
#endif  // CONTAINERS_EXPRESSIONS_ASDIAGONALMATRIX_CONTAINERS_ASDIAGONALMATRIX_TRAITS_HPP_
