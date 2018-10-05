
/*
//@HEADER
// ************************************************************************
//
//                                insert.hpp                               
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


#ifndef TINYMPL_INSERT_HPP
#define TINYMPL_INSERT_HPP

#include "variadic/erase.hpp"
#include "as_sequence.hpp"
#include "sequence.hpp"
#include "join.hpp"
#include "size.hpp"

namespace tinympl {

/**
 * \ingroup SeqAlgsIntr
 * \class insert
 * \brief Insert a subsequence into a given sequence at a given position
 * \param Pos The insertion position in the main sequence
 * \param SubSeq The subsequence
 * \param Seq The main sequence
 * \param Out The output sequence type
 */
template < std::size_t Pos,
  class SubSeq,
  class Seq,
  template<class ...> class Out = as_sequence<Seq>::template rebind>
struct insert : insert<Pos, as_sequence_t<SubSeq>, as_sequence_t<Seq>, Out> {};

template< std::size_t Pos, class ... SubSeqArgs, class ... SeqArgs,
  template<class...> class Out>
class insert<Pos, sequence<SubSeqArgs...>, sequence<SeqArgs...>, Out> {
    template<class ... Us>
    using head_seq = sequence<Us ..., SubSeqArgs ... >;

    typedef typename variadic::erase<Pos, sizeof ...( SeqArgs ), head_seq,
        SeqArgs ... >::type head;
    typedef typename variadic::erase<0, Pos, sequence, SeqArgs ... >::type tail;

  public:
    typedef typename join<Out<>, head, tail>::type type;
};

template <typename Sequence, typename Arg>
using push_front = insert<0, sequence<Arg>, Sequence>;
template <typename Sequence, typename Arg>
using push_front_t = typename push_front<Sequence, Arg>::type;

template <typename Sequence, typename SubSeq>
using extend_front = insert<0, as_sequence_t<SubSeq>, Sequence>;

template <typename Sequence, typename Arg>
using push_back = insert<size<Sequence>::value, sequence<Arg>, Sequence>;
template <typename Sequence, typename Arg>
using push_back_t = typename  push_back<Sequence, Arg>::type;

template <typename Sequence, typename SubSeq>
using extend_back = insert<size<Sequence>::value, as_sequence_t<SubSeq>, Sequence>;


} // namespace tinympl

#endif // TINYMPL_INSERT_HPP
