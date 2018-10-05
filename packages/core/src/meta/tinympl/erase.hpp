
/*
//@HEADER
// ************************************************************************
//
//                                erase.hpp                                
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


#ifndef TINYMPL_ERASE_HPP
#define TINYMPL_ERASE_HPP

#include "variadic/erase.hpp"
#include "as_sequence.hpp"
#include "sequence.hpp"
#include "transform.hpp"
#include "at.hpp"
#include "find_all_if.hpp"

namespace tinympl {

/**
 * \ingroup SeqAlgsIntr
 * \class erase
 * \brief Remove a range in a given sequence
 * \param First The first element to be removed
 * \param Last The first element which is not removed
 * \param Seq The input sequence
 * \param Out The output sequence type
 */
template<
  std::size_t First,
  std::size_t Last,
  class Seq,
  template <class...> class Out = as_sequence<Seq>::template rebind
>
struct erase : erase<First, Last, as_sequence_t<Seq>, Out> { };

template<
  std::size_t First,
  std::size_t Last,
  class... Args,
  template <class...> class Out
>
struct erase<First, Last, sequence<Args...>, Out> :
  variadic::erase<First, Last, Out, Args...> { };

template <
  class Seq,
  template <class...> class F,
  template <class...> class Out = as_sequence<Seq>::template rebind
>
struct erase_if_not {
  private:

    template <typename WrappedIndex>
    using _seq_at = at<WrappedIndex::value, Seq>;

  public:

    using type = tinympl::transform_t<
      tinympl::find_all_if_t<Seq, F>, _seq_at, Out
    >;
};

template <
  class Seq,
  template <class...> class F,
  template <class...> class Out = as_sequence<Seq>::template rebind
>
using erase_if_not_t = typename erase_if_not<Seq, F, Out>::type;

template <
  class Seq,
  template <class...> class F,
  template <class...> class Out = as_sequence<Seq>::template rebind
>
using erase_unless = erase_if_not<Seq, F, Out>;

template <
  class Seq,
  template <class...> class F,
  template <class...> class Out = as_sequence<Seq>::template rebind
>
using erase_unless_t = typename erase_if_not<Seq, F, Out>::type;

template <
  class Seq,
  template <class...> class F,
  template <class...> class Out = as_sequence<Seq>::template rebind
>
struct erase_if : erase_if_not<Seq, negate_metafunction<F>::template apply, Out> { };

template <
  class Seq,
  template <class...> class F,
  template <class...> class Out = as_sequence<Seq>::template rebind
>
using erase_if_t = typename erase_if<Seq, F, Out>::type;


} // namespace tinympl

#endif // TINYMPL_ERASE_HPP
