/*
//@HEADER
// ************************************************************************
//
// solvers_system_traits.hpp
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

#ifndef SOLVERS_EXPERIMENTAL_SYSTEM_TRAITS_HPP
#define SOLVERS_EXPERIMENTAL_SYSTEM_TRAITS_HPP

#include <type_traits>

#include "../../../containers/src/meta/containers_meta_detection_idiom.hpp"


namespace pressio{
namespace solvers{
namespace details {


template <typename T>
using has_public_matrix_type = typename T::matrix_type;


template <typename T>
using has_public_vector_type = typename T::vector_type;


template <
  typename T,
  typename Arg
>
using has_residual_callable_with_one_arg =
  decltype(std::declval<T>().residual(std::declval<Arg const&>()));


template <
  typename T,
  typename FirstArg,
  typename SecondArg
>
using has_residual_callable_with_two_args =
  decltype(std::declval<T>().residual(std::declval<FirstArg const&>(),
				      std::declval<SecondArg&>()));


template <
  typename T,
  typename Arg
>
using has_jacobian_callable_with_one_arg =
  decltype(std::declval<T>().jacobian(std::declval<Arg const&>()));


template <
  typename T,
  typename FirstArg,
  typename SecondArg
>
using has_jacobian_callable_with_two_args =
  decltype(std::declval<T>().jacobian(std::declval<FirstArg const&>(),
				      std::declval<SecondArg&>()));


template <typename T>
struct system_traits {
  typedef typename containers::meta::detected_t<has_public_vector_type, T> vector_type;
  typedef typename containers::meta::detected_t<has_public_matrix_type, T> matrix_type;

  static constexpr bool has_public_vector_type =
    containers::meta::is_detected<has_public_vector_type, T>::value;
  static constexpr bool has_public_matrix_type =
    containers::meta::is_detected<has_public_matrix_type, T>::value;

  static constexpr bool has_residual_callable_with_one_arg =
    containers::meta::is_detected<has_residual_callable_with_one_arg,
			    T, vector_type>::value;
  static constexpr bool has_residual_callable_with_two_args =
    containers::meta::is_detected<has_residual_callable_with_two_args,
			    T, vector_type, vector_type>::value;
  static constexpr bool has_residual_methods =
    has_residual_callable_with_one_arg || has_residual_callable_with_two_args;

  static constexpr bool has_jacobian_callable_with_one_arg =
    containers::meta::is_detected<has_jacobian_callable_with_one_arg,
			    T, vector_type>::value;
  static constexpr bool has_jacobian_callable_with_two_args =
    containers::meta::is_detected<has_jacobian_callable_with_two_args,
			    T, vector_type, matrix_type>::value;
  static constexpr bool has_jacobian_methods =
    has_jacobian_callable_with_one_arg || has_jacobian_callable_with_two_args;

  static constexpr bool is_system = has_residual_methods and  has_jacobian_methods;
};

} // end namespace details
} // end namespace solvers

}//end namespace pressio
#endif
