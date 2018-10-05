/*
//@HEADER
// ************************************************************************
//
//                     find_all.hpp
//                         darma
//              Copyright (C) 2015 Sandia Corporation
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


#ifndef TINYMPL_FIND_ALL_HPP
#define TINYMPL_FIND_ALL_HPP

#include "variadic/find_all.hpp"
#include "as_sequence.hpp"
#include "sequence.hpp"

namespace tinympl {

/**
 * \ingroup SeqNonModAlgs
 * \class find_all Compute the index of all elements in the sequence which are equal to the
 * type T.
 *
 * \param Sequence The input sequence
 * \param T The type to be tested
 * \return `find_all<...>::type` is `vector_c<std::size_t,v...>` where
`v...` are the 0-based indices of the elements equal_to T
 * \note The comparison is done with \ref tinympl::equal_to - it can be
specialized
 * \sa variadic::find_all
 */
template <class Sequence, class T>
struct find_all
  : public find_all<as_sequence_t<Sequence>, T>
{ };

template <class T, class... Args>
struct find_all<sequence<Args...>, T>
  : public variadic::find_all<T, Args...>
{ };

} // namespace tinympl

#endif // TINYMPL_FIND_ALL_HPP
