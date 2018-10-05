/*
//@HEADER
// ************************************************************************
//
//                          to_string.hpp
//                         darma_new
//              Copyright (C) 2017 NTESS, LLC
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

#ifndef SRC_META_TINYMPL_VARIADIC_TO_VECTOR_C_HPP_
#define SRC_META_TINYMPL_VARIADIC_TO_VECTOR_C_HPP_

#include "../string.hpp"
#include "../vector.hpp"
#include "all_of.hpp"
#include "../lambda.hpp"

namespace tinympl {
namespace variadic {

template <typename ValueType, typename... Args>
struct to_vector_c {
  private:
    typedef tinympl::vector<typename Args::value_type...> _value_type_vector;

    template <typename T>
    using _check = std::is_convertible<ValueType, typename T::value_type>;

    static_assert(all_of<_check, Args...>::value,
      "to_vector_c: unable to convert to string when all arguments do not have a value_type convertible to ValueType"
    );

  public:

    typedef vector_c<ValueType, (ValueType)Args::value...> type;
};

} // end namespace variadic
} // end namespace tinympl



#endif /* SRC_META_TINYMPL_VARIADIC_TO_VECTOR_C_HPP_ */
