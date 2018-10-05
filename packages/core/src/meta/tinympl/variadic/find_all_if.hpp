/*
//@HEADER
// ************************************************************************
//
//                          find_all_if.hpp
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

#ifndef SRC_META_TINYMPL_VARIADIC_FIND_ALL_IF_HPP_
#define SRC_META_TINYMPL_VARIADIC_FIND_ALL_IF_HPP_

#include <type_traits> // std::integral_constant, etc
#include <utility> // std::index_sequence

#include "../join.hpp"
#include "../transform.hpp"
#include "../as_sequence.hpp"
#include "../lambda.hpp"
#include "../bind.hpp"
#include "../plus.hpp"
#include "../to_index_sequence.hpp"
#include "../stl_integer_sequence.hpp"
#include "../logical_not.hpp"

namespace tinympl {
namespace variadic {

template <
  template <class...> class UnaryPredicate,
  typename... Args
>
struct find_all_if;

template <
  template <class...> class UnaryPredicate,
  typename Arg1, typename... Args
>
struct find_all_if<UnaryPredicate, Arg1, Args...> {
  private:
    template <typename T>
    using _add_one = typename plus<
      T,
      std::integral_constant<std::size_t, (std::size_t)1>
    >::type;

  public:
    using type = typename join<
      typename std::conditional_t<
        UnaryPredicate<Arg1>::type::value,
        std::index_sequence<0>,
        std::index_sequence<>
      >,
      typename tinympl::to_index_sequence<
        typename tinympl::transform<
          typename find_all_if<UnaryPredicate, Args...>::type,
          _add_one
        >::type // end transform
      >::type // end to_index_sequence
    >::type;
};

template <
  template <class...> class UnaryPredicate
>
struct find_all_if<UnaryPredicate> {
  using type = std::index_sequence<>;
};

template <
  template <class...> class UnaryPredicate,
  typename... Args
>
struct find_all_if_not
  : public find_all_if<
      negate_metafunction<UnaryPredicate>::template apply,
      Args...
    >
{ };


} // end namespace variadic
} // end namespace tinympl



#endif /* SRC_META_TINYMPL_VARIADIC_FIND_ALL_IF_HPP_ */
