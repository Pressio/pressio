/*
//@HEADER
// ************************************************************************
//
//                      partition.hpp
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

#ifndef TINYMPL_VARIADIC_PARTITION_HPP
#define TINYMPL_VARIADIC_PARTITION_HPP

#include <cstddef>

#include <tinympl/variadic/at.hpp>
#include <tinympl/splat.hpp>
#include <tinympl/join.hpp>
#include <tinympl/transform.hpp>
#include <tinympl/range_c.hpp>
#include <tinympl/sequence.hpp>
#include <tinympl/as_sequence.hpp>

namespace tinympl {
namespace variadic {

namespace _impl {

template <
  std::size_t NPerGroup, std::size_t IGroup, std::size_t NGroups,
  template <typename...> class InnerOut,
  template <typename...> class OuterOut,
  typename... Args
>
struct _partition {
  template <typename size_t_wrapped>
  using get_arg_at = at<size_t_wrapped::value, Args...>;

  using type = typename tinympl::join<
    OuterOut<
      typename tinympl::splat_to<
        typename tinympl::transform<
          tinympl::as_sequence_t<
            typename tinympl::make_range_c<
              std::size_t, IGroup*NPerGroup, (IGroup+1)*NPerGroup
            >::type
          >,
          get_arg_at
        >::type,
        InnerOut
      >::type
    >,
    typename _partition<
      NPerGroup, IGroup+1, NGroups, InnerOut, OuterOut, Args...
    >::type
  >::type;
};

template <
  std::size_t NPerGroup, std::size_t NGroups,
  template <typename...> class InnerOut,
  template <typename...> class OuterOut,
  typename... Args
>
struct _partition<NPerGroup, NGroups, NGroups, InnerOut, OuterOut, Args...> {
  using type = OuterOut<>;
};

} // end namespace _impl

template <
  std::size_t NPerGroup,
  template <typename...> class InnerOut,
  template <typename...> class OuterOut,
  typename... Args
>
struct partition
  : _impl::_partition<
      NPerGroup, 0ul, sizeof...(Args)/NPerGroup,
      InnerOut, OuterOut,
      Args...
    >
{
  static_assert(sizeof...(Args) % NPerGroup == 0,
    "invalid partition would create unequal sublists"
  );
};


} // end namespace variadic
} // end namespace tinympl


#endif //TINYMPL_VARIADIC_PARTITION_HPP
