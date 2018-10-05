
/*
//@HEADER
// ************************************************************************
//
//                              functional.hpp                             
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


#ifndef TINYMPL_FUNCTIONAL_HPP
#define TINYMPL_FUNCTIONAL_HPP

#include "int.hpp"
#include "short.hpp"
#include "char.hpp"
#include "bool.hpp"
#include "long.hpp"
#include "plus.hpp"
#include "multiplies.hpp"
#include "minus.hpp"
#include "divides.hpp"
#include "modulus.hpp"
#include "negate.hpp"
#include "equal_to.hpp"
#include "not_equal_to.hpp"
#include "less.hpp"
#include "greater.hpp"
#include "less_equal.hpp"
#include "greater_equal.hpp"
#include "and_b.hpp"
#include "or_b.hpp"
#include "not_b.hpp"
#include "logical_and.hpp"
#include "logical_or.hpp"
#include "logical_not.hpp"
#include "identity.hpp"
#include "inherit.hpp"
#include "sizeof.hpp"
#include "if.hpp"
#include "apply.hpp"

namespace tinympl {

/**
 * \defgroup NewTypes Type wrappers
 * Templates which wrap a value into a type.
 */

/**
 * \defgroup Functional Metafunctions
 * Class templates which implement metafunctions
 * @{
 */

  /**
   * \defgroup Arithmetic Arithmetic operations
   * Metafunctions which perform arithmetic operations on `std::integral_constant` or equivalent types
   */

  /**
   * \defgroup Comparisons Comparisons
   * Metafunctions which perform comparisons operations on `std::integral_constant` or equivalent types
   */

  /**
   * \defgroup Logical Logical operations
   * Metafunctions which perform logical operations on `std::integral_constant` or equivalent types
   */

/** @} */

}

#endif // TINYMPL_FUNCTIONAL_HPP
