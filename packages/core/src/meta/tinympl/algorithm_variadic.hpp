
/*
//@HEADER
// ************************************************************************
//
//                          algorithm_variadic.hpp                         
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


#ifndef TINYMPL_ALGORITHM_VARIADIC_HPP
#define TINYMPL_ALGORITHM_VARIADIC_HPP

#include "variadic/accumulate.hpp"
#include "variadic/all_of.hpp"
#include "variadic/any_of.hpp"
#include "variadic/copy.hpp"
#include "variadic/copy_if.hpp"
#include "variadic/copy_n.hpp"
#include "variadic/count.hpp"
#include "variadic/count_if.hpp"
#include "variadic/fill_n.hpp"
#include "variadic/find.hpp"
#include "variadic/find_if.hpp"
#include "variadic/generate_n.hpp"
#include "variadic/is_unique.hpp"
#include "variadic/left_fold.hpp"
#include "variadic/max_element.hpp"
#include "variadic/min_element.hpp"
#include "variadic/none_of.hpp"
#include "variadic/remove.hpp"
#include "variadic/remove_if.hpp"
#include "variadic/replace.hpp"
#include "variadic/replace_if.hpp"
#include "variadic/reverse.hpp"
#include "variadic/right_fold.hpp"
#include "variadic/sort.hpp"
#include "variadic/transform.hpp"
#include "variadic/unique.hpp"

// For backward compatibility
#include "bind.hpp"
#include "functional.hpp"
#include "variadic.hpp"

namespace tinympl {
namespace variadic {

/**
 * \defgroup VarAlgs Variadic algorithms
 * Algorithms which operate on variadic templates
 * @{
 */

  /**
   * \defgroup VarNonModAlgs Non-modifying sequence operations
   * Algorithms which analyze a sequence without producing an output sequence
   */

  /**
   * \defgroup VarModAlgs Modifying sequence operations
   * Algorithms which produce an output sequence
   */

  /**
   * \defgroup VarMaxMin Minimum/maximum operations
   * Algorithms which compute the minimum/maximum of a sequence
   */

  /**
   * \defgroup VarSort Sorting operations
   * Algorithms to sort a sequence.
   */

  /**
   * \defgroup VarSet Set operations (on unsorted sequences)
   * Algorithms which perform set operations.
   * \note Unlike the `std` counterparts, these algorithms do not require an
  ordering of the elements.
   */

/**
 * \defgroup VarFold Folding operations
 * Algorithms which perform reduction operations on a sequence.
 */

/** @} */


}
}

#endif // TINYMPL_ALGORITHM_VARIADIC_HPP
