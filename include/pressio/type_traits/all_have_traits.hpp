/*
//@HEADER
// ************************************************************************
//
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

#ifndef TYPE_TRAITS_ALL_HAVE_TRAITS_HPP_
#define TYPE_TRAITS_ALL_HAVE_TRAITS_HPP_

namespace pressio{ namespace impl{

template <class T, class = void>
struct has_traits : std::false_type{};

template <class T>
struct has_traits<
  T, std::enable_if_t<
       // I need to check for something meaninful
       // because if I only check for ::pressio::Traits<T>
       // any type would yield true since ::pressio::Traits<T>
       // is always the default case
       mpl::not_void<
	 typename ::pressio::Traits<T>::scalar_type
	 >::value
       &&
       ::pressio::Traits<T>::rank != 0
       >
  > : std::true_type{};

} //end namespace impl

template <class ... Args>
struct all_have_traits;

template <class T>
struct all_have_traits<T>{
  static constexpr auto value = impl::has_traits<T>::value;
};

template <class T1, class T2>
struct all_have_traits<T1, T2>{
  static constexpr auto value = impl::has_traits<T1>::value
    && impl::has_traits<T2>::value;
};

template <class T1, class T2, class T3, class ... rest>
struct all_have_traits<T1, T2, T3, rest...>{
  static constexpr auto value =
    all_have_traits<T1, T2>::value &&
    all_have_traits<T3, rest...>::value;
};

}
#endif  // TYPE_TRAITS_ALL_HAVE_TRAITS_HPP_
