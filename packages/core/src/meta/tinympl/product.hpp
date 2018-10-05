/*
//@HEADER
// ************************************************************************
//
//                      product.hpp
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

#ifndef DARMAFRONTEND_PRODUCT_HPP
#define DARMAFRONTEND_PRODUCT_HPP

#include "as_sequence.hpp"
#include "vector.hpp"

namespace tinympl {

namespace _impl {

template <
  template <class...> class F,
  typename NewArgLists,
  typename OldArgLists,
  typename... RemainingSeqs
>
struct product_helper;

template <
  template <class...> class F,
  typename NewArgLists,
  typename OldArgLists,
  typename Sequence,
  typename... RemainingSequences
>
struct product_helper<
  F, NewArgLists, OldArgLists, Sequence, RemainingSequences...
> : product_helper<F, NewArgLists, OldArgLists, tinympl::as_sequence_t<Sequence>, RemainingSequences...>
{ };


template <
  template <class...> class F,
  typename NewArgLists,
  typename SequenceEntry,
  typename... OldArgListsEntries,
  typename... RemainingSequenceEntries,
  typename... RemainingSequences
>
struct product_helper<
  F, NewArgLists,
  vector<OldArgListsEntries...>,
  sequence<SequenceEntry, RemainingSequenceEntries...>,
  RemainingSequences...
>
{
  using type = typename product_helper<
    F,
    // For each of the old entries, extend the new list to include a list
    // with SequenceEntry
    typename NewArgLists::template extend<vector<
      typename OldArgListsEntries::template push_back<SequenceEntry>::type...
    >>::type,
    // Pass through the old list
    vector<OldArgListsEntries...>,
    // and pop off a sequence element
    sequence<RemainingSequenceEntries...>,
    RemainingSequences...
  >::type;
};

template <
  template <class...> class F,
  typename NewArgLists,
  typename OldArgLists,
  typename... RemainingSequences
>
struct product_helper<
  F, NewArgLists,
  OldArgLists,
  sequence<>,
  RemainingSequences...
>
{
  using type = typename product_helper<
    F,
    // Reset the new arg lists
    tinympl::vector<>,
    // move the new arg lists to the old arg lists
    NewArgLists,
    // and pop off the empty sequence
    RemainingSequences...
  >::type;
};

template <
  template <class...> class F,
  typename OldArgLists,
  typename... FinalArgListsEntries
>
struct product_helper<
  F, vector<FinalArgListsEntries...>,
  OldArgLists,
  sequence<>
>
{
  using type = tinympl::vector<
    splat_to_t<FinalArgListsEntries, F>...
  >;
};


} // end namespace _impl

template <
  template <class...> class F,
  class... Sequences
>
struct product
  : _impl::product_helper<F, vector<>, vector<vector<>>, Sequences...>
{ };

} // end namespace tinympl

#endif //DARMAFRONTEND_PRODUCT_HPP
