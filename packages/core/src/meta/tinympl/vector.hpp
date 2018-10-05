
/*
//@HEADER
// ************************************************************************
//
//                                vector.hpp                               
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


#ifndef TINYMPL_VECTOR_HPP
#define TINYMPL_VECTOR_HPP

#include "variadic/at.hpp"
#include "variadic/erase.hpp"
#include "erase.hpp"
#include "insert.hpp"
#include "splat.hpp"
#include "find_if.hpp"


namespace tinympl {

/**
 * \defgroup Containers Containers
 * Full and half compile time containers of types and values.
 * @{
 */

/**
 * \class vector
 * \brief A compile time vector of types
 * Vector is the simplest tinympl sequence type. It provides standard modifiers and random access
 * to its elements.
 */
template <class... Args>
struct vector
{
  enum {
    size = sizeof...(Args) //!< The size of the vector
  };

  enum {
    empty = (size == 0) //!< Determine whether the vector is empty
  };

  //! Access the i-th element
  template <std::size_t i>
  struct at {
    static_assert(i < size, "Index i is out of range");
    using type = typename variadic::at<i,Args...>::type;
  };
  template <std::size_t i>
  using at_t = typename at<i>::type;

  //! Access the i-th element, but safe protected
  template <typename Default, std::size_t i>
  struct at_or {
    using type = typename variadic::at_or<Default, i, Args...>::type;
  };
  template<typename Default, std::size_t i>
  using at_or_t = typename at_or<Default, i>::type;

  //! Return a new vector constructed by inserting `T` on the back of the current vector
  template <class T>
  struct push_back {
    using type = vector<Args..., T>;
  };
  template <class T>
  using push_back_t = typename push_back<T>::type;

  //! Return a new vector constructed by inserting `T` on the front of the current vector
  template <class T>
  struct push_front {
    using type = vector<T, Args...>;
  };
  template <class T>
  using push_front_t = typename push_front<T>::type;

  //! Return a new vector constructed by removing the last element of the current vector
  struct pop_back {
    using type = typename variadic::erase<size-1,size,tinympl::vector,Args...>::type;
  };
  struct safe_pop_back {
    using type = typename std::conditional_t<
      (size > 0),
      variadic::erase<size-1, size, tinympl::vector, Args...>,
      identity<vector<>>
    >::type;
  };
  // Can't/Shouldn't have a pop_back_t because it won't be lazily constructed

  //! Return a new vector constructed by removing the first element of the current vector
  struct pop_front {
    using type = typename variadic::erase<0, 1, tinympl::vector,Args...>::type;
  };
  // Can't/Shouldn't have a pop_front_t because it won't be lazily constructed

  //! Return a new vector constructed by erasing the elements in the range [first,last) of the current vector
  template <std::size_t first, std::size_t last>
  struct erase : tinympl::erase<first, last, vector<Args...>, tinympl::vector> {};
  template <std::size_t first, std::size_t last>
  using erase_t = typename erase<first, last>::type;

  template <template <class...> class F>
  struct erase_if : tinympl::erase_if<vector<Args...>, F, tinympl::vector> {};
  template <template <class...> class F>
  using erase_if_t = typename erase_if<F>::type;

  template <template <class...> class F>
  struct erase_if_not : tinympl::erase_if_not<vector<Args...>, F, tinympl::vector> {};
  template <template <class...> class F>
  using erase_if_not_t = typename erase_if_not<F>::type;
  template <template <class...> class F>
  using erase_unless = erase_if_not<F>;
  template <template <class...> class F>
  using erase_unless_t = typename erase_if_not<F>::type;

  template <typename Seq>
  struct extend : tinympl::join<vector<Args...>, Seq> { };
  template <typename Seq>
  using extend_t = typename extend<Seq>::type;
  template <typename Seq>
  using extend_back = extend<Seq>;
  template <typename Seq>
  using extend_back_t = typename extend<Seq>::type;

  template <typename Seq>
  struct extend_front
    : tinympl::join<
        tinympl::splat_to_t<Seq, tinympl::vector>,
        tinympl::vector<Args...>
      >
  { };
  template <typename Seq>
  using extend_front_t = typename extend_front<Seq>::type;

  template <template <class...> class F>
  using splat_to = tinympl::splat_to<sequence<Args...>, F>;
  template <template <class...> class F>
  using splat_to_t = tinympl::splat_to_t<sequence<Args...>, F>;

  //! Return a new vector constructed by inserting the elements `Ts...` in the current vector starting at the index `i`
  template <std::size_t i, class... Ts>
  struct insert : tinympl::insert<i,
    sequence<Ts...>,
    vector<Args...>,
    tinympl::vector
  > { };
  template <std::size_t i, class... Ts>
  using insert_t = typename insert<i, Ts...>::type;

  template <template <class...> class F, typename Default>
  struct get_first_if_or_default
    : tinympl::at_or<
        Default,
        tinympl::find_if_t<sequence<Args...>, F>::value,
        sequence<Args...>
      >
  { };

  template <template <class...> class F, typename Default>
  struct get_last_if_or_default
    : tinympl::at_or<
        Default,
        tinympl::find_last_if_t<sequence<Args...>, F>::value,
        sequence<Args...>
      >
  { };

  //! Return the first element of the vector
  struct front { using type = variadic::at_t<0, Args...>; };
  // Can't/Shouldn't have a front_t because it won't be lazily constructed

  //! Return the last element of the vector
  struct back { using type = variadic::at_t<size-1, Args...>; };
  // Can't/Shouldn't have a back_t because it won't be lazily constructed

  //! Return the last element of the vector
  template <typename Default = tinympl::nonesuch>
  struct safe_back {
    using type = typename std::conditional_t<
      (size > 0),
      variadic::at<size-1, Args...>,
      identity<Default>
    >::type;
  };
};

/** @} */

/**
 * \ingroup SeqCustom
 * \brief Customization point to allow `vector` to work as a tinympl sequence
 */
template <class... Args>
struct as_sequence<vector<Args...> > {
  using type = sequence<Args...>;
  template <class... Ts> using rebind = vector<Ts...>;
};

}

#endif // TINYMPL_VECTOR_HPP
