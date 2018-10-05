/*
//@HEADER
// ************************************************************************
//
//                      switch.hpp
//                         DARMA
//              Copyright (C) 2017 Sandia Corporation
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

#ifndef TINYMPL_SWITCH_HPP
#define TINYMPL_SWITCH_HPP

#include <tinympl/equal_to.hpp>

namespace tinympl {

// Workaround for ICC issue with default template template parameters
namespace _impl {


} // end namespace _impl

template <
  typename ValueType,
  ValueType value,
  typename DefaultResult = void,
  // metafunction class with apply taking two ValueType arguments
  typename value_equal_mfc = tinympl::make_value_equal<ValueType>
>
struct static_value_switch {
  template <
    ValueType case_value,
    typename Result,
    typename FoundResult,
    bool result_found = false
  >
  struct _value_switch_case {
    using type = std::conditional_t<
      result_found,
      FoundResult,
      std::conditional_t<
        value_equal_mfc::template apply<value, case_value>::type::value,
        Result,
        DefaultResult
      >
    >;
    template <
      ValueType next_value,
      typename NextResult
    >
    using case_ = _value_switch_case<
      next_value, NextResult,
      type,
      result_found || value_equal_mfc::template apply<value, case_value>::type::value
    >;
    template <
      typename DefResult
    >
    using default_ = _value_switch_case<
      case_value, /* ignored */
      Result, /* ignored */
      std::conditional_t<result_found, FoundResult, DefResult>,
      true
    >;
  };
  template <
    ValueType next_value,
    typename NextResult
  >
  using case_ = _value_switch_case<
    next_value, NextResult,
    DefaultResult, false
  >;

};

template <
  typename ValueType,
  ValueType value,
  typename DefaultResult = void,
  typename value_equal = tinympl::make_value_equal<ValueType>
>
using switch_value_ = static_value_switch<
  ValueType, value, DefaultResult, value_equal
>;

template <
  typename MainType,
  typename DefaultResult,
  template <class...> class Equal = tinympl::equal_to
>
struct static_switch {
  template <
    typename CompareTo,
    typename Result,
    typename FoundResult,
    bool result_found = false
  >
  struct _value_switch_case {
    using type = std::conditional_t<
      result_found,
      FoundResult,
      std::conditional_t<
        Equal<MainType, CompareTo>::type::value,
        Result,
        DefaultResult
      >
    >;
    template <
      typename NextCompareTo,
      typename NextResult
    >
    using case_ = _value_switch_case<
      NextCompareTo, NextResult, type,
      result_found || Equal<MainType, CompareTo>::type::value
    >;
    template <
      typename DefResult
    >
    using default_ = _value_switch_case<
      CompareTo, /* ignored */
      Result, /* ignored */
      std::conditional_t<result_found, FoundResult, DefResult>,
      true
    >;
  };
  template <
    typename CompareTo,
    typename NextResult
  >
  using case_ = _value_switch_case<
    CompareTo, NextResult,
    DefaultResult, false
  >;

};

} // end namespace tinympl


#endif //TINYMPL_SWITCH_HPP
