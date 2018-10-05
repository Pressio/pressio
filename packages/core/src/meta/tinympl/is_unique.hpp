
/*
//@HEADER
// ************************************************************************
//
//                              is_unique.hpp                              
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


#ifndef TINYMPL_IS_UNIQUE_HPP
#define TINYMPL_IS_UNIQUE_HPP

#include "sequence.hpp"
#include "as_sequence.hpp"
#include "variadic/is_unique.hpp"

namespace tinympl {

/**
 * \ingroup SeqSet
 * \class is_unique
 * \brief Determines whether the input sequence contains only unique elements
 * \param Sequence the input sequence
 * \return `is_unique<...>::type` is a `std::integral_constant<bool,v>` where
`v` is true iff the input sequence contains no duplicates
 * \note Unlike `std::sort`, the input sequence is not required to be sorted,
but the compile time complexity is O(N^2)
 * \note The comparison is done with \ref tinympl::equal_to - it can be
specialized
 * \sa variadic::is_unique
 */
template<class Sequence> struct is_unique :
    is_unique<as_sequence_t<Sequence> > {};

template<class ... Args> struct is_unique<sequence<Args...> > :
    variadic::is_unique<Args...> {};

} // namespace tinympl

#endif // TINYMPL_IS_UNIQUE_HPP
