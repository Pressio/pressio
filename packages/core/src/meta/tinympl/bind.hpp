
/*
//@HEADER
// ************************************************************************
//
//                                 bind.hpp                                
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


#ifndef TINYMPL_BIND_HPP
#define TINYMPL_BIND_HPP

#include "variadic.hpp"

namespace tinympl {

/**
 * \defgroup BindAndLambda Bind and Lambda expressions.
 * Helpers to transform and define metafunction classes on the fly.
 * @{
 */

/**
 * \class bind
 * \brief Produce a template type by binding the given arguments on the passed template template.
 *
 * `tinympl::bind` is the compile time equivalent of `std::bind`. It produces a new template
 * type `bind<...>::template eval` which invokes the given one (`F`) with some of its arguments bound to `Args`.
 * Notice that in C++11 the effect of bind can be achieved with template aliases. In order to produce a cleaner
 * code, we recommend to use template aliases wherever is possible, and use `bind` only when necessary.
 *
 * `bind< std::is_same, arg1, U>::template eval` is equivalent to
 * `template<class T> using is_same1 = std::is_same<T,U>;`
 *
 * `bind` also automatically recognize bind expressions in its subarguments, so it is possible to nest multiple bind calls:
 *
 * `bind< std::is_same, bind<std::remove_reference,arg1>, U>::template eval` is equivalent to
 * `template<class T> using is_same1 = std::is_same< typename std::remove_reference<T>::type, U>;`
 */
template< template<class ... T> class F,class ... Args> struct bind;

template<std::size_t> struct arg;
typedef arg<1> arg1;
typedef arg<2> arg2;
typedef arg<3> arg3;
typedef arg<4> arg4;
typedef arg<5> arg5;
typedef arg<6> arg6;
typedef arg<7> arg7;
typedef arg<8> arg8;

namespace placeholders {

typedef arg<1> _1;
typedef arg<2> _2;
typedef arg<3> _3;
typedef arg<4> _4;
typedef arg<5> _5;
typedef arg<6> _6;
typedef arg<7> _7;
typedef arg<8> _8;
typedef arg<9> _9;
typedef arg<1> _;

} // end namespace placeholders

/**
 * \brief Determine whether a type is a placeholder.
 * `is_placeholder<T>::value` is 0 if `T` is not a placeholder, otherwise is the index of the placeholder
 */
template <class T>
struct is_placeholder
  : std::integral_constant<std::size_t, 0> { };

template<std::size_t i>
struct is_placeholder< arg<i> >
  : std::integral_constant<std::size_t, i>
{
  static_assert(i != 0, "Placeholder arg<0> is undefined");
  template <typename T> struct cv_qualifier_rebind { typedef T type; };
  //template <typename T> using cv_qualifier_rebind_t =
  //  typename cv_qualifier_rebind<T>::type;
};

template<std::size_t i>
struct is_placeholder< const arg<i> >
  : std::integral_constant<std::size_t, i>
{
  static_assert(i != 0, "Placeholder arg<0> is undefined");
  template <typename T> struct cv_qualifier_rebind { typedef const T type; };
  //template <typename T> using cv_qualifier_rebind_t =
  //  typename cv_qualifier_rebind<T>::type;
};

template<std::size_t i>
struct is_placeholder< const volatile arg<i> >
  : std::integral_constant<std::size_t, i>
{
  static_assert(i != 0, "Placeholder arg<0> is undefined");
  template <typename T> struct cv_qualifier_rebind { typedef const volatile T type; };
  //template <typename T> using cv_qualifier_rebind_t =
  //  typename cv_qualifier_rebind<T>::type;
};

template<std::size_t i>
struct is_placeholder< volatile arg<i> >
  : std::integral_constant<std::size_t, i>
{
  static_assert(i != 0, "Placeholder arg<0> is undefined");
  template <typename T> struct cv_qualifier_rebind { typedef volatile T type; };
  //template <typename T> using cv_qualifier_rebind_t =
  //  typename cv_qualifier_rebind<T>::type;
};

/**
 * \brief Determine whether a type is a bind expression.
 */
template<class T> struct is_bind_expression : std::false_type {};
template<template<class ... T> class F,class ... Args> struct is_bind_expression< bind<F,Args...> > : std::true_type {};

/** @} */

template <template <class... T> class F, class Head, class... Tail>
struct bind<F, Head, Tail...>
{
  private:
    template <class... Args>
    struct call
    {
      template <class... BoundArgs>
      struct eval
      {
        template<class T,class Enable = void> struct pick {typedef T type;};
        template<class T> struct pick<T, typename std::enable_if< (is_placeholder<T>::value > 0) >::type> {typedef variadic::at_t<is_placeholder<T>::value-1, Args ... > type;};
        template<class T> struct pick<T, typename std::enable_if< is_bind_expression<T>::value >::type> {typedef typename T::template eval<Args...>::type type;};

        typedef typename pick<Head>::type argument_t;

        //Forward the call to bind
        typedef typename bind<F,Tail...>::template call<Args...>::template eval<BoundArgs..., argument_t>::type type;
      };
    };

    template< template<class ...> class,class ...> friend struct bind;

  public:
    template<class ... Args>
    struct eval
    {
      using type = typename call<Args...>::template eval<>::type;
    };

    template<class... Ts>
    struct eval_value
    {
      static constexpr auto value = eval<Ts...>::type::value;
    };

    template<class ... Args>
    using eval_t = typename eval<Args...>::type;
    template<class ... Args>
    using apply = eval<Args...>;
    template<class ... Args>
    using apply_value = eval_value<Args...>;
};

template <template <class... T> class F>
struct bind<F>
{
  private:
    template <class... Args>
    struct call
    {
      template <class... BoundArgs>
      struct eval
      {
          typedef typename F<BoundArgs...>::type type;
      };
    };

    template <template <class...> class, class...> friend struct bind;

  public:
    template <class... Args>
    struct eval
    {
      using type = typename call<Args...>::template eval<>::type;
    };
    template<class... Ts>
    struct eval_value
    {
      static constexpr auto value = eval<Ts...>::type::value;
    };

    template<class ... Args>
    using eval_t = typename eval<Args...>::type;
    template<class ... Args>
    using apply = eval<Args...>;
    template<class ... Args>
    using apply_value = eval_value<Args...>;
};

}

#endif // MPL_BIND_HPP
