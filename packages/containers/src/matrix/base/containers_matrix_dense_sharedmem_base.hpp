/*
//@HEADER
// ************************************************************************
//
// containers_matrix_dense_sharedmem_base.hpp
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

#ifndef CONTAINERS_MATRIX_BASE_MATRIX_DENSE_SHAREDMEM_BASE_HPP_
#define CONTAINERS_MATRIX_BASE_MATRIX_DENSE_SHAREDMEM_BASE_HPP_

#include "../containers_matrix_traits.hpp"

namespace pressio{
namespace containers{

template<typename derived_type>
class MatrixDenseSharedMemBase
{
  static_assert( details::traits<derived_type>::is_shared_mem==1,
  "OOPS: NON-shared mem matrix inheriting from dense sharedMem base!");

  using this_t = MatrixDenseSharedMemBase<derived_type>;
  using traits = ::pressio::containers::details::traits<derived_type>;
  using subspan_ret_t = typename traits::subspan_ret_t;
  using subspan_const_ret_t = typename traits::subspan_const_ret_t;

public:
  template< typename _subspan_ret_t = subspan_ret_t>
  mpl::enable_if_t< !std::is_void<_subspan_ret_t>::value, _subspan_ret_t>
  subspan(const typename _subspan_ret_t::interval_t & rowRangeIn,
	  const typename _subspan_ret_t::interval_t & colRangeIn)
  {
    auto & matObj = static_cast<derived_type &>(*this);
    return subspan_ret_t(matObj, rowRangeIn, colRangeIn);
  };

  template< typename _subspan_const_ret_t = subspan_const_ret_t>
  mpl::enable_if_t< !std::is_void<_subspan_const_ret_t>::value, _subspan_const_ret_t>
  subspan(const typename _subspan_const_ret_t::interval_t & rowRangeIn,
	  const typename _subspan_const_ret_t::interval_t & colRangeIn) const
  {
    const auto & matObj = static_cast<const derived_type &>(*this);
    return subspan_const_ret_t(matObj, rowRangeIn, colRangeIn);
  };

private:
  friend derived_type;
  using sc_t = typename details::traits<derived_type>::scalar_t;
  using ord_t = typename details::traits<derived_type>::ordinal_t;

};//end class

} // end namespace containers
}//end namespace pressio
#endif
