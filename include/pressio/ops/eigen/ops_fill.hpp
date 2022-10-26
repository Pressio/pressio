/*
//@HEADER
// ************************************************************************
//
// ops_fill.hpp
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

#ifndef OPS_EIGEN_OPS_FILL_HPP_
#define OPS_EIGEN_OPS_FILL_HPP_

namespace pressio{ namespace ops{

/*constrained via is_convertible because the impl is using
  native Eigen ops which use expressions, so value has to
  be convertible to the the scalar type of the object */
template <typename T, class ScalarType>
#ifdef PRESSIO_ENABLE_CXX20
requires (::pressio::is_native_container_eigen<T>::value
       || ::pressio::is_expression_acting_on_eigen<T>::value)
       && std::convertible_to<ScalarType, typename ::pressio::Traits<T>::scalar_type>
void
#else
::pressio::mpl::enable_if_t<
    (::pressio::is_native_container_eigen<T>::value
  || ::pressio::is_expression_acting_on_eigen<T>::value)
  && std::is_convertible<ScalarType, typename ::pressio::Traits<T>::scalar_type>::value
  >
#endif
fill(T & o, const ScalarType & value)
{
  const typename ::pressio::Traits<T>::scalar_type v(value);
  impl::get_native(o).setConstant(v);
}

}}//end namespace pressio::ops
#endif  // OPS_EIGEN_OPS_FILL_HPP_
