/*
//@HEADER
// ************************************************************************
//
//                          copy_traits.hpp
//                         darma_new
//              Copyright (C) 2017 NTESS, LLC
//
// Under the terms of Contract DE-NA-0003525 with NTESS, LLC,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact somebody@sandia.gov
//
// ************************************************************************
//@HEADER
*/

#ifndef SRC_META_TINYMPL_COPY_TRAITS_HPP_
#define SRC_META_TINYMPL_COPY_TRAITS_HPP_

namespace tinympl {

#define _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES \
  template <typename T> using apply_t = typename apply<T>::type; \
  template <typename T> using rebind = apply<T>; \
  template <typename T> using rebind_t = typename apply<T>::type;


//template <typename From>
//struct copy_cv_qualifiers<const From> {
//  template <typename T> struct apply { typedef const T type; };
//  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
//};
//
//template <typename From>
//struct copy_cv_qualifiers<const volatile From> {
//  template <typename T> struct apply { typedef const volatile T type; };
//  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
//};
//
//template <typename From>
//struct copy_cv_qualifiers<volatile From> {
//  template <typename T> struct apply { typedef volatile T type; };
//  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
//};

template <typename From>
struct copy_volatileness {
  template <typename T> struct apply { typedef T type; };
  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
};

template <typename From>
struct copy_volatileness<volatile From> {
  template <typename T> struct apply { typedef volatile T type; };
  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
};

template <typename From>
struct copy_constness {
  template <typename T> struct apply { typedef T type; };
  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
};

template <typename From>
struct copy_constness<const From> {
  template <typename T> struct apply { typedef const T type; };
  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
};

template <typename From>
struct copy_cv_qualifiers {
  template <typename T> struct apply {
    typedef typename copy_constness<From>::template apply<
      typename copy_volatileness<From>::template apply<T>::type
    >::type type;
  };
  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
};

template <typename From>
struct copy_extents {
  template <typename T> struct apply { typedef T type; };
  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
};

template <typename From>
struct copy_extents<From[]> {
  template <typename T> struct apply { typedef T type[]; };
  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
};

template <typename From, size_t N>
struct copy_extents<From[N]> {
  template <typename T> struct apply { typedef T type[N]; };
  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
};

template <typename From>
struct copy_all_extents {
  template <typename T> struct apply { typedef T type; };
  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
};

template <typename From>
struct copy_all_extents<From[]> {
  template <typename T> struct apply {
    typedef typename copy_all_extents<From>::template apply<T>::type type[];
  };
  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
};

template <typename From, size_t N>
struct copy_all_extents<From[N]> {
  template <typename T> struct apply {
    typedef typename copy_all_extents<From>::template apply<T>::type type[N];
  };
  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
};

template <typename From>
struct copy_reference_type {
  template <typename T> struct apply { typedef typename std::remove_reference<T>::type type; };
  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
};

template <typename From>
struct copy_reference_type<From&> {
  template <typename T> struct apply { typedef typename std::remove_reference<T>::type& type; };
  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
};

template <typename From>
struct copy_reference_type<From&&> {
  template <typename T> struct apply { typedef typename std::remove_reference<T>::type&& type; };
  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
};

template <typename From>
struct copy_all_type_properties {
  template <typename T> struct apply {
    typedef typename copy_cv_qualifiers<From>::template apply<
      typename copy_all_extents<From>::template apply<
        typename copy_reference_type<From>::template apply<T>::type
      >::type
    >::type type;
  };
  _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES
};


#undef _TINYMPL_tmp_COPY_TRAITS_APPLY_ALIASES

} // end namespace tinympl

#endif /* SRC_META_TINYMPL_COPY_TRAITS_HPP_ */
