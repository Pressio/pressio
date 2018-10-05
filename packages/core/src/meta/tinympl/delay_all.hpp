/*
//@HEADER
// ************************************************************************
//
//                      delay_all.hpp
//                         DARMA
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

#ifndef TINYMPL_DELAY_ALL_HPP
#define TINYMPL_DELAY_ALL_HPP

#include <tinympl/delay.hpp>
#include <tinympl/vector.hpp>

namespace tinympl {

template <typename... Args>
struct delayed_arguments {
  template <template <class...> class Metafunction>
  struct apply_metafunction {
    using type = typename Metafunction<Args...>::type;
  };
  template <typename MetafunctionClass>
  struct apply_metafunction_class {
    using type = typename MetafunctionClass::template apply<Args...>::type;
  };
  using type = tinympl::vector<Args...>;
};



namespace _impl {

template <typename T>
using _is_delayed_application_archetype = bool_< T::_tinympl__is_delayed_application >;

template <typename T>
using is_delayed_application = detected_or<std::false_type,
  _is_delayed_application_archetype, T
>;

template <typename OuterMFC, typename UnwrappedArgsVector, typename... Args>
struct _delay_all_helper;

template <
  typename OuterMFC, typename UnwrappedArgsVector,
  typename InnerMFCOrNormalArg, typename PotentiallyWrappedDelayedArgs,
  typename... Args
>
struct _delay_all_helper<
  OuterMFC, UnwrappedArgsVector,
  InnerMFCOrNormalArg, PotentiallyWrappedDelayedArgs, Args...
>
{
  template <typename NormalArg, typename NextArg, typename Enable=void>
  struct next_helper {
    using type = typename _delay_all_helper<OuterMFC,
      push_back_t<UnwrappedArgsVector, NormalArg>, NextArg, Args...
    >::type;
  };

  template <typename DelayedArg, typename NextArg>
  struct next_helper<DelayedArg, NextArg,
    std::enable_if_t<is_delayed_application<DelayedArg>::value>
  >
  {
    using type = typename _delay_all_helper<OuterMFC,
      push_back_t<UnwrappedArgsVector, typename DelayedArg::type>,
    NextArg, Args...
    >::type;
  };

  template <typename InnerMFC, typename... WrappedArgs>
  struct next_helper<InnerMFC, delayed_arguments<WrappedArgs...>> {
    // TODO decide if it's a good idea to static_assert here, since it short-circuits SFINAE
    static_assert(
      is_metafunction_class<InnerMFC>::value,
      "delayed_arguments given without a valid metafunction class preceding it"
    );
    using type = typename _delay_all_helper<OuterMFC,
      push_back_t<UnwrappedArgsVector,
      typename _delay_all_helper<InnerMFC, tinympl::vector<>,
        WrappedArgs...
      >::type
    >, Args...
    >::type;
  };

  using type = typename next_helper<
    InnerMFCOrNormalArg, PotentiallyWrappedDelayedArgs
  >::type;

};

template < typename OuterMFC, typename UnwrappedArgsVector, typename LastArg >
struct _delay_all_helper<
  OuterMFC, UnwrappedArgsVector, LastArg
>
{
  template <typename NormalArg, typename Enable=void>
  struct next_helper {
    using type = typename _delay_all_helper<OuterMFC,
      push_back_t<UnwrappedArgsVector, NormalArg>
    >::type;
  };

  template <typename DelayedArg>
  struct next_helper<DelayedArg,
    std::enable_if_t<is_delayed_application<DelayedArg>::value>
  >
  {
    using type = typename _delay_all_helper<OuterMFC,
      push_back_t<UnwrappedArgsVector, typename DelayedArg::type>
    >::type;
  };

  using type = typename next_helper<LastArg>::type;
};

template <
  typename OuterMFC, typename... UnwrappedArgs
>
struct _delay_all_helper<OuterMFC, tinympl::vector<UnwrappedArgs...>> {
  using type = typename OuterMFC::template apply<UnwrappedArgs...>::type;
};


} // end namespace _impl

template <typename MetafunctionClass, typename... WrappedArgs>
struct delay_all {
  using type = typename _impl::_delay_all_helper<
    MetafunctionClass, tinympl::vector<>,
    WrappedArgs...
  >::type;
};

template <template <class...> class F>
struct delayed_metafunction {
  template <typename... Args>
  struct applied_to {
    static constexpr auto _tinympl__is_delayed_application = true;
    using type = typename delay_all<
      metafunction_wrap<F>, Args...
    >::type;
  };
  template <typename... Args> using on = applied_to<Args...>;
  template <typename... Args>
  struct instantiated_with {
    static constexpr auto _tinympl__is_delayed_application = true;
    using type = typename delay_all<
      wrap_instantiation<F>, Args...
    >::type;
  };
  template <typename... Args> using of = instantiated_with<Args...>;
};
template <template <class...> class F>
using lazy = delayed_metafunction<F>;

template <typename T>
struct delayed {
  struct applied {
    using type = typename T::type;
  };
  using lazy_get_type = applied;
};


template <
  template <class...> class F,
  typename... Types
>
struct delay_instantiate {
  using type = F<
    typename Types::type...
  >;
};

template <
  template <class...> class Delayed
>
struct undelay {
  template <typename... Args>
  struct apply {
    typedef typename Delayed<Args...>::type::type type;
  };
};

namespace _impl {

template <typename T>
using _is_lazy_value_suspension_archetype = decltype( T::type::value );

template <typename T>
using _has_member_value_archetype = decltype( T::value );

} // end namespace _impl

template <typename T>
using is_lazy_value_suspension = is_detected<
  _impl::_is_lazy_value_suspension_archetype, T>;

template <typename T>
using is_lazy_boolean_value_suspension = is_detected_convertible<bool,
  _impl::_is_lazy_value_suspension_archetype, T>;

template <typename T>
using has_member_value = is_detected<_impl::_has_member_value_archetype, T>;

template <typename T>
using has_member_value_convertible_to_bool =
is_detected_convertible<bool, _impl::_has_member_value_archetype, T>;

template <typename T>
struct extract_value_potentially_lazy {
  private:
    using _type_to_evaluate =
    typename std::conditional_t<
      has_member_value<T>::value,
      identity<T>, T
    >::type;

  public:

    static constexpr auto value = _type_to_evaluate::value;

    using type = std::integral_constant<decltype(value), value>;

};
template <typename T>
using lazy_get_value = extract_value_potentially_lazy<T>;

template <typename T>
struct extract_bool_value_potentially_lazy {
  private:
    using _type_to_evaluate = typename std::conditional_t<
      has_member_value_convertible_to_bool<T>::value,
      identity<T>, T
    >::type;

  public:

    static constexpr bool value = _type_to_evaluate::value;

    using type = std::integral_constant<bool, value>;
};

template <typename T>
struct extract_value_potentially_lazy_arbitrary_depth
{
  private:

    template <typename U>
    struct _get_type_to_eval {
      using type = typename std::conditional_t<
        has_member_value<U>::value,
        identity<U>, _get_type_to_eval<typename U::type>
      >::type;
    };

    using _type_to_evaluate = typename _get_type_to_eval<T>::type;

  public:

    static constexpr auto value = _type_to_evaluate::value;

    using type = std::integral_constant<decltype(value), value>;

};

template <typename T>
struct extract_bool_value_potentially_lazy_arbitrary_depth
{
  private:

    template <typename U>
    struct _get_type_to_eval {
      using type = typename std::conditional_t<
        has_member_value_convertible_to_bool<U>::value,
        identity<U>, _get_type_to_eval<typename U::type>
      >::type;
    };

    using _type_to_evaluate = typename _get_type_to_eval<T>::type;

  public:

    static constexpr bool value = _type_to_evaluate::value;

    using type = std::integral_constant<bool, value>;

};


} // end namespace tinympl

#endif // TINYMPL_DELAY_ALL_HPP
