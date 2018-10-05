
/*
//@HEADER
// ************************************************************************
//
//                              algorithm.hpp                              
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


#ifndef TINYMPL_ALGORITHM_HPP
#define TINYMPL_ALGORITHM_HPP

#include "accumulate.hpp"
#include "all_of.hpp"
#include "any_of.hpp"
#include "at.hpp"
#include "copy.hpp"
#include "copy_if.hpp"
#include "copy_n.hpp"
#include "count.hpp"
#include "count_if.hpp"
#include "erase.hpp"
#include "fill_n.hpp"
#include "find.hpp"
#include "find_if.hpp"
#include "generate_n.hpp"
#include "insert.hpp"
#include "is_unique.hpp"
#include "join.hpp"
#include "left_fold.hpp"
#include "lexicographical_compare.hpp"
#include "max_element.hpp"
#include "min_element.hpp"
#include "none_of.hpp"
#include "remove.hpp"
#include "remove_if.hpp"
#include "replace.hpp"
#include "replace_if.hpp"
#include "reverse.hpp"
#include "right_fold.hpp"
#include "set_difference.hpp"
#include "set_intersection.hpp"
#include "set_union.hpp"
#include "size.hpp"
#include "sort.hpp"
#include "transform.hpp"
#include "transform2.hpp"
#include "transform_many.hpp"
#include "transpose.hpp"
#include "unique.hpp"
#include "unordered_equal.hpp"
#include "zip.hpp"

// For backward compatibility.
#include "algorithm_variadic.hpp"
#include "functional.hpp"
#include "sequence.hpp"

namespace tinympl {

/**
 * \defgroup SeqAlgs Sequence algorithms
 * Algorithms which operate on sequences
 * @{
 */

  /**
   * \defgroup SeqAlgsIntr Sequence introspection operations
   * Algorithms which perform trivial operations on sequences
   */

  /**
   * \defgroup SeqNonModAlgs Non-modifying sequence operations
   * Algorithms which analyze a sequence without producing an output sequence
   */

  /**
   * \defgroup SeqModAlgs Modifying sequence operations
   * Algorithms which produce an output sequence
   */

  /**
   * \defgroup SeqMaxMin Minimum/maximum operations
   * Algorithms which compute the minimum/maximum of a sequence and
  lexicographical compare
   */

  /**
   * \defgroup SeqSort Sorting operations
   * Algorithms to sort a sequence.
   */

  /**
   * \defgroup SeqSet Set operations (on unsorted sequences)
   * Algorithms which perform set operations.
   * \note Unlike the `std` counterparts, these algorithms do not require an
  ordering of the elements.
   */

  /**
   * \defgroup SeqFold Folding operations
   * Algorithms which perform reduction operations on a sequence.
   */

/** @} */

}

#endif // TINYMPL_ALGORITHM_HPP
