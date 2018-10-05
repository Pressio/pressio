
/*
//@HEADER
// ************************************************************************
//
//                              generate_n.hpp                             
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


#ifndef TINYMPL_GENERATE_N_HPP
#define TINYMPL_GENERATE_N_HPP

#include "variadic/generate_n.hpp"

namespace tinympl {

/**
 * \ingroup SeqModAlgs
 * \class generate_n
 * \brief Generate N elements using a given generator metafunction
 * \param N The number of elements to generate
 * \param Gen The generator. `Gen< std::integral_constant<int,i> >::type` must
be a valid expression.
 * \param Out the output sequence type
 * \return `generate_n<...>::type` is a type templated from `Out` constructed
with N elements generated with `Gen< int_<0> >, Gen< int_<1> >, ... Gen<
int_<N-1> >`
 * \sa variadic::generate_n
 */
template<std::size_t N,
        template<class ...> class Gen,
        template<class ...> class Out>
struct generate_n : variadic::generate_n<N, Gen, Out> {};

} // namespace tinympl

#endif // TINYMPL_GENERATE_N_HPP
