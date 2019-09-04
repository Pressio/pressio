/*
//@HEADER
// ************************************************************************
//
// solvers_meta_static_checks.hpp
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

#ifndef SOLVERS_META_META_STATIC_CHECKS_HPP
#define SOLVERS_META_META_STATIC_CHECKS_HPP

#include <type_traits>

#include "../../../containers/src/matrix/containers_matrix_traits.hpp"
#include "../../../containers/src/vector/containers_vector_traits.hpp"


namespace pressio{
namespace solvers{
namespace meta {

/**
 * @brief Check whether two matrices are compatible.
 *
 * @section DESCRIPTION
 *
 * Two matrices are compatible iff the following two conditions are satisfied:
 *   1. the underlying matrix representation is the same; this implies that the
 *      underlying linear containers library is also the same;
 *   2. the matrices have the same structure (sparse or dense).
 */
template <typename T, typename U>
struct are_matrix_compatible {
  static constexpr bool valid_matrix = containers::details::traits<T>::wrapped_package_identifier != containers::details::WrappedPackageIdentifier::Undefined;
  static constexpr bool value = valid_matrix && (containers::details::traits<T>::wrapped_matrix_identifier == containers::details::traits<U>::wrapped_matrix_identifier);
};


/**
 * @brief Check whether a vector and a matrix are compatible.
 *
 * @section DESCRIPTION
 *
 * A vector and a matrix are compatible iff the underlying linear containers package
 * used to represent them is the same.
 */
template <typename T, typename U>
struct are_vector_matrix_compatible {
  static constexpr bool valid_vector = containers::details::traits<T>::wrapped_package_identifier != containers::details::WrappedPackageIdentifier::Undefined;
  static constexpr bool valid_matrix = containers::details::traits<U>::wrapped_package_identifier != containers::details::WrappedPackageIdentifier::Undefined;
  static constexpr bool value = valid_vector && valid_matrix && (containers::details::traits<T>::wrapped_package_identifier == containers::details::traits<U>::wrapped_package_identifier);
};


/**
 * @brief Check whether two vectors are compatible.
 *
 * @section DESCRIPTION
 *
 * Two vectors are compatible iff their underlying structure is the same and they
 * are represented using the same underlying linear containers package.
 */
template <
  typename T,
  typename U,
  typename Enable = void
>
struct are_vector_compatible;


/**
 * Either one or both arguments are not vectors.
 */
template <
  typename T,
  typename U
>
struct are_vector_compatible<
  T,
  U,
  typename std::enable_if<
    !containers::details::traits<T>::is_vector || !containers::details::traits<U>::is_vector,
    void
  >::type
> {
  static constexpr bool value = false;
};


/**
 * Both arguments are vectors.
 */
template <
  typename T,
  typename U
>
struct are_vector_compatible<
  T,
  U,
  typename std::enable_if<
    containers::details::traits<T>::is_vector && containers::details::traits<U>::is_vector,
    void
  >::type
> {

  static constexpr bool valid_vector =
    containers::details::traits<T>::wrapped_package_identifier !=
    containers::details::WrappedPackageIdentifier::Undefined;

  static constexpr bool same_type =
    valid_vector &&
    containers::details::traits<T>::wrapped_vector_identifier ==
    containers::details::traits<U>::wrapped_vector_identifier;

  static constexpr bool value = same_type;
    // (containers::details::traits<T>::is_dynamic ||
    //  containers::details::traits<U>::is_dynamic ||
    //  containers::details::traits<T>::rows == containers::details::traits<U>::rows);
};


} // end namespace meta
} // end namespace solvers
}//end namespace pressio
#endif
