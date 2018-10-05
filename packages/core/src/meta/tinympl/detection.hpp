/*
//@HEADER
// ************************************************************************
//
//                          detection.hpp
//                         tinympl
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

#ifndef TINYMPL_DETECTION_H_
#define TINYMPL_DETECTION_H_

#include <type_traits>

#include "void_t.hpp"

namespace tinympl {

// Large pieces taken or adapted from http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2015/n4436.pdf

// primary template handles all types not supporting the archetypal Op
template <
  class Default,
  class _always_void,
  template <class...> class Op,
  class... Args
>
struct detector {
  constexpr static auto value = false;
  using type = Default;
};

// specialization recognizes and handles only types supporting Op
template <
  class Default,
  template <class...> class Op,
  class... Args
>
struct detector<Default, void_t<Op<Args...>>, Op, Args...> {
  constexpr static auto value = true;
  using type = Op<Args...>;
};

struct nonesuch {
  nonesuch() = delete;
  ~nonesuch() = delete;
  nonesuch(nonesuch const&) = delete;
  void operator=(nonesuch const&) = delete;
};

// TODO make this a tinympl::bool_
template <template <class...> class Op, class... Args>
using is_detected = detector<nonesuch, void, Op, Args...>;

template <template <class...> class Op, class... Args>
using detected_t = typename is_detected<Op, Args...>::type;

template <class Default, template <class...> class Op, class... Args>
using detected_or = detector<Default, void, Op, Args...>;

template <class Default, template <class...> class Op, class... Args>
using detected_or_t = typename detected_or<Default, Op, Args...>::type;

template <class Expected, template<class...> class Op, class... Args>
using is_detected_exact = std::is_same<Expected, detected_t<Op, Args...>>;

template <class To, template <class...> class Op, class... Args>
using is_detected_convertible = std::is_convertible<detected_t<Op, Args...>, To>;

} // end namespace tinympl

#define _TINYMPL_DETECTED_TYPE_IMPL( \
  assign_to, detected_name, of_type, default_type, access_specifier \
) \
  private: \
    template <typename _##detected_name##_arch_param> \
    using _##assign_to##_detect_##detected_name##_archetype = \
      typename _##detected_name##_arch_param::detected_name; \
  access_specifier: \
    using assign_to = ::tinympl::detected_or_t<default_type, \
       _##assign_to##_detect_##detected_name##_archetype, \
       of_type \
    >

#define TINYMPL_PUBLIC_DETECTED_TYPE(assign_to, detected_name, of_type) \
  _TINYMPL_DETECTED_TYPE_IMPL( \
    assign_to, detected_name, of_type, ::tinympl::nonesuch, public \
  )

#define TINYMPL_PUBLIC_DETECTED_TYPE_WITH_DEFAULT( \
  assign_to, detected_name, of_type, default_type \
) \
  _TINYMPL_DETECTED_TYPE_IMPL( \
    assign_to, detected_name, of_type, default_type, public \
  )

#define TINYMPL_PRIVATE_DETECTED_TYPE(assign_to, detected_name, of_type) \
  _TINYMPL_DETECTED_TYPE_IMPL( \
    assign_to, detected_name, of_type, ::tinympl::nonesuch, private \
  )

#define TINYMPL_PRIVATE_DETECTED_TYPE_WITH_DEFAULT( \
  assign_to, detected_name, of_type, default_type \
) \
  _TINYMPL_DETECTED_TYPE_IMPL( \
    assign_to, detected_name, of_type, default_type, private \
  )

#define TINYMPL_PROTECTED_DETECTED_TYPE(assign_to, detected_name, of_type) \
  _TINYMPL_DETECTED_TYPE_IMPL( \
    assign_to, detected_name, of_type, ::tinympl::nonesuch, protected \
  )

#define TINYMPL_PROTECTED_DETECTED_TYPE_WITH_DEFAULT( \
assign_to, detected_name, of_type, default_type \
) \
  _TINYMPL_DETECTED_TYPE_IMPL( \
    assign_to, detected_name, of_type, default_type, protected \
  )

#define _TINYMPL_DETECTED_VALUE_WITH_DEFAULT_IMPL( \
  assign_to, detected_name, of_type, default_value, access_specifier \
) \
  private: \
    template <typename _##detected_name##_arch_param> \
    using _##assign_to##_detect_value_##detected_name##_archetype = \
      ::std::integral_constant< \
        decltype(_##detected_name##_arch_param::detected_name), \
        _##detected_name##_arch_param::detected_name \
      >; \
  access_specifier: \
    static constexpr auto assign_to = ::tinympl::detected_or_t< \
       ::std::integral_constant<decltype(default_value), default_value>, \
       _##assign_to##_detect_value_##detected_name##_archetype, \
       of_type \
    >::value

#define TINYMPL_PUBLIC_DETECTED_VALUE(assign_to, detected_name, of_type) \
  _TINYMPL_DETECTED_VALUE_WITH_DEFAULT_IMPL( \
    assign_to, detected_name, of_type, ::std::false_type, public \
  )

#define TINYMPL_PUBLIC_DETECTED_VALUE_WITH_DEFAULT( \
  assign_to, detected_name, of_type, default_value \
) \
  _TINYMPL_DETECTED_VALUE_WITH_DEFAULT_IMPL( \
    assign_to, detected_name, of_type, default_value, public \
  )

#define TINYMPL_PRIVATE_DETECTED_VALUE(assign_to, detected_name, of_type) \
  _TINYMPL_DETECTED_VALUE_WITH_DEFAULT_IMPL( \
    assign_to, detected_name, of_type, ::std::false_type, private \
  )

#define TINYMPL_PRIVATE_DETECTED_VALUE_WITH_DEFAULT( \
assign_to, detected_name, of_type, default_value \
) \
  _TINYMPL_DETECTED_VALUE_WITH_DEFAULT_IMPL( \
    assign_to, detected_name, of_type, default_value, private \
  )

#define TINYMPL_PROTECTED_DETECTED_VALUE(assign_to, detected_name, of_type) \
  _TINYMPL_DETECTED_VALUE_WITH_DEFAULT_IMPL( \
    assign_to, detected_name, of_type, ::std::false_type, protected \
  )

#define TINYMPL_PROTECTED_DETECTED_VALUE_WITH_DEFAULT( \
assign_to, detected_name, of_type, default_value \
) \
  _TINYMPL_DETECTED_VALUE_WITH_DEFAULT_IMPL( \
    assign_to, detected_name, of_type, default_value, protected \
  )

#endif /* TINYMPL_DETECTION_H_ */
