/*
//@HEADER
// ************************************************************************
//
//                          is_suspension.hpp
//                         whatever
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

#ifndef META_TINYMPL_IS_SUSPENSION_HPP_
#define META_TINYMPL_IS_SUSPENSION_HPP_

#include <type_traits>

#include "is_instantiation.hpp"
#include "logical_and.hpp"
#include "detection.hpp"

namespace tinympl {

namespace _impl {

template <typename T>
using _has_member_named_type_archetype = typename T::type;

template <typename T>
using _has_member_type_named_type =
  detected_or<nonesuch, _has_member_named_type_archetype, T>;

} // end namespace _impl


template <typename T>
struct is_suspension
  : and_<
      is_instantiation<T>,
      _impl::_has_member_type_named_type<T>
    >
{ };

namespace _impl {

template <typename T>
using has_value_archetype = decltype( T::value );

template <typename T>
using has_value = is_detected_convertible<bool, has_value_archetype, T>;

} // end namespace _impl


template <typename T>
struct is_value_suspension
  : and_<
    is_instantiation<T>,
    _impl::has_value<T>
  >
{ };

} // end namespace tinympl



#endif /* META_TINYMPL_IS_SUSPENSION_HPP_ */
