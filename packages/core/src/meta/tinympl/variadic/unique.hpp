/*
//@HEADER
// ************************************************************************
//
//                                unique.hpp                               
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


#ifndef TINYMPL_VARIADIC_UNIQUE_HPP
#define TINYMPL_VARIADIC_UNIQUE_HPP

#include "remove.hpp"

namespace tinympl {
namespace variadic {

/**
 * \ingroup VarModAlgs
 * \class unique
 * \brief Produces a sequence of unique elements from the input sequence,
preserving the ordering.
 * \param Out the output sequence type - defaults to the same kind of the input
sequence.
 * \param Args... The input sequence.
 * \return `unique<...>::type` is a type templated from `Out` which contains the
resulting sequence.
 * \note Only the first (leftmost) duplicate is mantained in the output
sequence.
 * \sa tinympl::unique
 */
template<template<class ...> class Out, class ... Args> struct unique;

template<template<class ...> class Out, class Head, class ... Tail>
struct unique<Out, Head, Tail...> {
private:
    template<class ... Ts>
    struct impl {
        template<class ... Us> using next = unique<Out, Us...>;
        typedef typename remove<Head, next, Tail...>::type::
            template impl<Ts..., Head>::type type;
    };

    template<template<class...> class, class...> friend struct unique;

public:
    typedef typename impl<>::type type;

};

template<template<class ...> class Out> struct unique<Out> {
private:
    template<class ... Ts>
    struct impl {
        typedef Out<Ts...> type;
    };

    template<template<class...> class, class...> friend struct unique;

public:
    typedef typename impl<>::type type;
};


} // namespace variadic
} // namespace tinympl

#endif // TINYMPL_VARIADIC_UNIQUE_HPP
