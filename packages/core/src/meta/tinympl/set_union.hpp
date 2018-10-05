
/*
//@HEADER
// ************************************************************************
//
//                              set_union.hpp                              
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


#ifndef TINYMPL_SET_UNION_HPP
#define TINYMPL_SET_UNION_HPP

#include "variadic/unique.hpp"
#include "as_sequence.hpp"
#include "sequence.hpp"

namespace tinympl {

/**
 * \ingroup SeqSet
 * \class set_union
 * \brief Computes the union of two sets
 * \param SequenceA The sequence which represents the first set - duplicates are
ignored
 * \param SequenceB The sequence which represents the second set - duplicates
are ignored
 * \param Out the output sequence type - defaults to the same kind of
`SequenceA`
 * \return A type templated from `Out` which contains the resulting sequence
 * \note If `seq<A...>` and `seq<B...>` are the two input sequences, the
resulting sequence is constructed as `unique< seq<A...,B...> >`
 */
template<class SequenceA,
        class SequenceB,
        template<class ...> class Out = as_sequence<SequenceA>::template rebind>
struct set_union :
    set_union< as_sequence_t<SequenceA>, as_sequence_t<SequenceB>, Out> {};

template<class ... Ts, class ... Us, template<class ...> class Out>
struct set_union<sequence<Ts...>, sequence<Us...>, Out> :
    variadic::unique<Out, Ts..., Us...> {};

} // namespace tinympl

#endif // TINYMPL_SET_UNION_HPP
