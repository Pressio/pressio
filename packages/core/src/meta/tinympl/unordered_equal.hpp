
/*
//@HEADER
// ************************************************************************
//
//                           unordered_equal.hpp                           
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


#ifndef TINYMPL_UNORDERED_EQUAL_HPP
#define TINYMPL_UNORDERED_EQUAL_HPP

#include "variadic/count.hpp"
#include "variadic/all_of.hpp"
#include "as_sequence.hpp"
#include "sequence.hpp"
#include "logical_and.hpp"
#include "variadic/count.hpp"
#include "variadic/all_of.hpp"
#include <type_traits>

namespace tinympl {

/**
 * \ingroup SeqSet
 * \class unordered_equal
 * \brief Determines whether it is possible to reorder the sequence `A` to match
exactly the sequence `B`
 * \param SequenceA The first sequence
 * \param SequenceB The second sequence
 * \return `unordered_equal<A,B>::type` is a `std::integral_constant<bool,v>`
where `v` is true only if the two sequences are equal (except for the ordering
of the elements)
 * \note Compile time complexity is O(N^2)
 */
template<class SequenceA, class SequenceB>
struct unordered_equal :
    unordered_equal< as_sequence_t<SequenceA>, as_sequence_t<SequenceB> > {};

namespace detail {
template<class SequenceA, class SequenceB> struct unordered_equal_impl;

template<class ... As, class ... Bs>
struct unordered_equal_impl<sequence<As...>, sequence<Bs...> > {
private:
    template<class T>
    struct check_t {
        typedef std::integral_constant < bool,
                variadic::count<T, As...>::type::value ==
                variadic::count<T, Bs...>::type::value > type;
    };

public:
    typedef typename logical_and <
        typename variadic::all_of< check_t, As...>::type,
        typename variadic::all_of< check_t, Bs...>::type >::type type;
};
}

template<class ... As, class ... Bs>
struct unordered_equal<sequence<As...>, sequence<Bs...> > :
    detail::unordered_equal_impl< sequence<As...>, sequence<Bs...> >::type
{};

} // namespace tinympl

#endif // TINYMPL_UNORDERED_EQUAL_HPP
