/*
//@HEADER
// ************************************************************************
//
// ops_matching_extents.hpp
//                     		  Pressio
//                             Copyright 2019
//    National Technology & Engineering Solutions of Sandia, LLC (NTESS)
//
// Under the terms of Contract DE-NA0003525 with NTESS, the
// U.S. Government retains certain rights in this software.
//
// Pressio is licensed under BSD-3-Clause terms of use:
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived
// from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
// STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
// IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Francesco Rizzi (fnrizzi@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef OPS_OPS_MATCHING_EXTENTS_HPP_
#define OPS_OPS_MATCHING_EXTENTS_HPP_

namespace pressio{ namespace ops{

template<class T1, class T2>
struct matching_extents<T1, T2>
{
  static bool compare(const T1 & a, const T2 & b)
  {
    using size_type = typename ::pressio::Traits<T1>::ordinal_type;
    for (size_type i=0; i<::pressio::Traits<T1>::rank; ++i)
      {
	if(::pressio::ops::extent(a, i) != ::pressio::ops::extent(b, i)){
	  return false;
	}
      }
    return true;
  }
};

template<class T1, class T2, class T3>
struct matching_extents<T1, T2, T3>
{
  static bool compare(const T1 & a, const T2 & b, const T3& c)
  {
    using size_type = typename ::pressio::Traits<T1>::ordinal_type;
    for (size_type i=0; i<::pressio::Traits<T1>::rank; ++i)
      {
	if(::pressio::ops::extent(a, i)!=::pressio::ops::extent(b, i) or
	   ::pressio::ops::extent(a, i)!=::pressio::ops::extent(c, i))
	  {
	    return false;
	  }
      }
    return true;
  }
};

template<class T1, class T2, class T3, class T4>
struct matching_extents<T1, T2, T3, T4>
{
  static bool compare(const T1 & a, const T2 & b, const T3& c, const T4& d)
  {
    using size_type = typename ::pressio::Traits<T1>::ordinal_type;
    for (size_type i=0; i<::pressio::Traits<T1>::rank; ++i)
      {
	if(::pressio::ops::extent(a, i)!=::pressio::ops::extent(b, i) or
	   ::pressio::ops::extent(a, i)!=::pressio::ops::extent(c, i) or
	   ::pressio::ops::extent(a, i)!=::pressio::ops::extent(d, i))
	  {
	    return false;
	  }
      }
    return true;
  }
};

template<class T1, class T2, class T3, class T4, class T5>
struct matching_extents<T1, T2, T3, T4, T5>
{
  static bool compare(const T1 & a, const T2 & b,
                      const T3& c, const T4& d,
                      const T5 & e)
  {
    using size_type = typename ::pressio::Traits<T1>::ordinal_type;
    for (size_type i=0; i<::pressio::Traits<T1>::rank; ++i){
      if(::pressio::ops::extent(a, i)!=::pressio::ops::extent(b, i) or
         ::pressio::ops::extent(a, i)!=::pressio::ops::extent(c, i) or
         ::pressio::ops::extent(a, i)!=::pressio::ops::extent(d, i) or
         ::pressio::ops::extent(a, i)!=::pressio::ops::extent(e, i))
        {
	  return false;
        }
    }
    return true;
  }
};

}}//end namespace pressio::ops
#endif  // OPS_OPS_MATCHING_EXTENTS_HPP_
