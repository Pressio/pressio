/*
//@HEADER
// ************************************************************************
//
//                     find_all_if.hpp
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


#ifndef TINYMPL_FIND_ALL_IF_HPP
#define TINYMPL_FIND_ALL_IF_HPP

#include "variadic/find_all_if.hpp"
#include "as_sequence.hpp"
#include "sequence.hpp"

namespace tinympl {

/**
 * \ingroup SeqNonModAlgs
 * \class find_all_if Compute the index of all elements in the sequence which
satisfies a given predicate, as a vector_c<std::size_t>
 * \param Sequence The input sequence
 * \param F The test predicate - `F<T>::type::value` shall be convertible to bool
 * \return `find_all_if<...>::type` is `vector_c<std::size_t,v...>` where
`v...` are the 0-based indices of the elements which satisfy `F`.
 * \sa variadic::find_all_if
 */
template <class Sequence, template <class... T> class F>
struct find_all_if
  : public find_all_if<as_sequence_t<Sequence>, F>
{ };

template <template <class... T> class F, class... Args>
struct find_all_if<sequence<Args...>, F>
  : public variadic::find_all_if<F, Args...>
{ };

template <class Sequence, template <class... T> class F>
using find_all_if_t = typename find_all_if<Sequence, F>::type;

} // namespace tinympl

#endif // TINYMPL_FIND_ALL_IF_HPP
