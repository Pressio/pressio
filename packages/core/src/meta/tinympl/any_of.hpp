
/*
//@HEADER
// ************************************************************************
//
//                                any_of.hpp                               
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


#ifndef TINYMPL_ANY_OF_HPP
#define TINYMPL_ANY_OF_HPP

#include <tinympl/variadic/any_of.hpp>
#include <tinympl/as_sequence.hpp>
#include <tinympl/sequence.hpp>

namespace tinympl {

/**
 * \ingroup SeqNonModAlgs
 * \class any_of
 * \brief Determines whether any of the elements in the sequence satisfy the
given predicate
 * \param Sequence the input sequence
 * \param F the predicate, `F<T>::type::value` must be convertible to `bool`
 * \return `any_of<...>::type` is a `std::integral_constant<bool,v>` where `v`
is true iff at least one element in the sequence satisfy the predicate `F`
 * \sa variadic::any_of
 */
template <class Sequence, template <class...> class F>
struct any_of : any_of<as_sequence_t<Sequence>, F> { };

template< template<class ...> class F, class ... Args>
struct any_of<sequence<Args...>, F > : variadic::any_of<F, Args...> { };

namespace types_only {

template< typename S, class FWrapped>
struct any_of : tinympl::any_of<
  S, FWrapped::template apply_value
> { };

} // end namespace types_only

} // namespace tinympl

#endif // TINYMPL_ANY_OF_HPP
