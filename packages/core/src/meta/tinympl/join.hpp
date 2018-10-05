
/*
//@HEADER
// ************************************************************************
//
//                                 join.hpp                                
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


#ifndef TINYMPL_JOIN_HPP
#define TINYMPL_JOIN_HPP

#include "sequence.hpp"
#include "as_sequence.hpp"

namespace tinympl {

/**
 * \ingroup SeqAlgsIntr
 * \class join
 * \brief Merge two sequences
 * \param Args The sequences
 * \return A sequence constructed by joining all the passed sequences with the
type of the first one
 */
template <class... Args> struct join;

template <class Head, class Next, class ... Tail>
struct join<Head, Next, Tail...> {
  typedef
    typename join<typename join<Head, Next>::type, Tail...>::type type;
};

template <class Head, class Last>
struct join<Head, Last> {
  private:
    template <class S1, class S2, template <class...> class Out>
    struct do_join;

    template <class... S1, class... S2, template <class...> class Out>
    struct do_join<sequence<S1...>, sequence<S2...>, Out> {
        typedef Out<S1..., S2...> type;
    };

  public:
    typedef typename do_join<
      as_sequence_t<Head>,
      as_sequence_t<Last>,
      as_sequence<Head>::template rebind
    >::type type;
};

template <class Head>
struct join<Head> {
  typedef Head type;
};

} // namespace tinympl

#endif // TINYMPL_JOIN_HPP
