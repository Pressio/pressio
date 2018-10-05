
/*
//@HEADER
// ************************************************************************
//
//                               identity.hpp                              
//                         whatever
//              Copyright (C) 2015 Sandia Corporation
// This file was adapted from its original form in the tinympl library.
// The original file bore the following copyright:
//   Copyright (C) 2013, Ennio Barbaro.
// See LEGAL.md for more information.
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


#ifndef TINYMPL_IDENTITY_HPP
#define TINYMPL_IDENTITY_HPP

namespace tinympl {

/**
 * \ingroup Functional
 * \class identity
 * \brief Returns the argument passed
 */
template<class T> struct identity {typedef T type;};

/**
 * \ingroup Functional
 * \class identity_value
 * \brief Returns the argument passed and T::value for value
 */
template<class T>
struct value_identity {
  typedef T type;
  static constexpr decltype(type::value) value = type::value;
};

/**
 * \ingroup Functional
 * \class ignore_argument
 * \brief ignores the first argument and instead returns the second argument,
 * which defaults to void
 */
template <typename /* ignored */, typename ReturnType=void>
struct ignore_argument {
  using type = ReturnType;
};

template <typename ReturnType>
struct make_ignore_argument {
  template <typename /* ignored */>
  using apply = tinympl::identity<ReturnType>;
};

/**
 * \ingroup Functional
 * \class ignore_value_argument
 * \brief ignores the non-type argument (of type T) and instead returns
 * the third argument, which defaults to void
 */
template <typename T, T /*ignored*/, typename ReturnType=void>
struct ignore_value_argument {
  using type = ReturnType;
};

} // namespace tinympl

#endif // TINYMPL_IDENTITY_HPP
