
/*
//@HEADER
// ************************************************************************
//
//                       lexicographical_compare.hpp                       
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


#ifndef TINYMPL_LEXICOGRAPHICAL_COMPARE_HPP
#define TINYMPL_LEXICOGRAPHICAL_COMPARE_HPP

#include "as_sequence.hpp"
#include "sequence.hpp"
#include "less.hpp"
#include <type_traits>

namespace tinympl {

/**
 * \ingroup SeqMaxMin
 * \class lexicographical_compare
 * \brief Compares two sequences using the given comparator
 * \param SequenceA the First sequence
 * \param SequenceB the second sequence
 * \param Comparator the comparator (default less)
 * \return An `std::integral_constant<int,v>` where `v` is -1 if the first
sequence is lexicographically smaller than the second, 1 if the first is greater
than the second and 0 if the two sequences are equal.
 */
template<class SequenceA,
        class SequenceB,
        template<class ...> class Comparator = less>
struct lexicographical_compare :
    lexicographical_compare<as_sequence_t<SequenceA>,
                            as_sequence_t<SequenceB>,
                            Comparator> {};

template<class AHead,
        class ... ATail,
        class BHead,
        class ... BTail,
        template<class ...> class Comparator>
struct lexicographical_compare<
    sequence<AHead, ATail...>,
    sequence<BHead, BTail...>,
    Comparator> :
        std::integral_constant < int,
        ( Comparator<AHead, BHead>::type::value ?
              -1 :
              Comparator<BHead, AHead>::type::value ?
              1 :
              lexicographical_compare <
                  sequence<ATail...>,
                  sequence<BTail...>,
                  Comparator >::value
        ) > {};

template<class ... As, template<class ...> class Comparator>
struct lexicographical_compare< sequence<As...>, sequence<>, Comparator> :
        std::integral_constant<int, 1> {};

template<class ... Bs, template<class ...> class Comparator>
struct lexicographical_compare< sequence<>, sequence<Bs...>, Comparator> :
        std::integral_constant < int, -1 > {};

template<template<class ...> class Comparator>
struct lexicographical_compare< sequence<>, sequence<>, Comparator> :
        std::integral_constant<int, 0> {};

} // namespace tinympl

#endif // TINYMPL_LEXICOGRAPHICAL_COMPARE_HPP
