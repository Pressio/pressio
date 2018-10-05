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

#ifndef TINYMPL_PARTITION_HPP
#define TINYMPL_PARTITION_HPP

#include <tinympl/variadic/partition.hpp>
#include <tinympl/as_sequence.hpp>
#include <tinympl/sequence.hpp>

namespace tinympl {

/**
 *  @brief e.g.,
 *    partition<2, Seq<Arg1, Arg2, Arg3, Arg4, Arg5, Arg6>,
 *      Inner, Outer
 *    >::type = Outer<Inner<Arg1, Arg2>, Inner<Arg3, Arg4>, Inner<Arg5, Arg6>>
 */
template <
  size_t NPerGroup,
  typename Sequence,
  template <class...> class InnerOut = as_sequence<Sequence>::template rebind,
  template <class...> class OuterOut = as_sequence<Sequence>::template rebind
>
struct partition
  : partition<NPerGroup, as_sequence_t<Sequence>, InnerOut, OuterOut> { };

template <
  size_t NPerGroup,
  template <typename...> class InnerOut,
  template <typename...> class OuterOut,
  typename... Args
>
struct partition<NPerGroup, sequence<Args...>, InnerOut, OuterOut>
  : tinympl::variadic::partition<NPerGroup, InnerOut, OuterOut, Args...> { };


} // end namespace tinympl

#endif //DARMA_PARTITION_HPP
