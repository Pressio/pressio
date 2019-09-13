/*
//@HEADER
// ************************************************************************
//
// rom_operator_base.hpp
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

#ifndef ROM_OPERATOR_BASE_HPP_
#define ROM_OPERATOR_BASE_HPP_

#include "rom_ConfigDefs.hpp"

namespace pressio{ namespace rom{

template <typename derived_t>
class OperatorBase
  : private utils::details::CrtpBase<OperatorBase<derived_t>>{

public:

  // A (whatever that is) acts on X, return result
  template <typename operand_t>
  auto apply(const operand_t & X)
    -> decltype(
    std::declval<derived_t>().template applyImpl<operand_t>(X) ){
    return this->underlying().applyImpl(X);
  }

  // A (whatever that is) acts on X
  // result stored in Y
  template <typename operand_t1,
	    typename operand_t2>
  void apply(const operand_t1 & X, operand_t2 & Y){
    this->underlying().applyImpl(X,Y);
  }
  //---------------------------------------------------------
  //---------------------------------------------------------

  // A (whatever that is) acts on X from right, return result
  // => X A
  template <typename operand_t>
  auto applyRight(const operand_t & X)
    -> decltype(
    std::declval<derived_t>().template applyRightImpl<operand_t>(X) ){
    return this->underlying().applyRightImpl(X);
  }

  // A (whatever that is) acts on X from right, return result
  // => X A
  template <typename operand_t1,
	    typename operand_t2>
  void applyRight(const operand_t1 & X, operand_t2 & Y){
    this->underlying().applyRightImpl(X, Y);
  }
  //---------------------------------------------------------
  //---------------------------------------------------------

  // A^T (whatever that is) acts on X, return result
  template <typename operand_t>
  auto applyTranspose(const operand_t & X)
    -> decltype(
    std::declval<derived_t>().template applyTransposeImpl<operand_t>(X) ){
    return this->underlying().applyTransposeImpl(X);
  }

  // A^T (whatever that is) acts on X,
  // result stored in Y
  template <typename operand_t1,
	    typename operand_t2>
  void applyTranspose(const operand_t1 & X,
		      operand_t2 & Y){
    this->underlying().applyTransposeImpl(X,Y);
  }
  //---------------------------------------------------------
  //---------------------------------------------------------

private:
  friend derived_t;
  friend utils::details::CrtpBase<OperatorBase<derived_t>>;
  OperatorBase() = default;
  ~OperatorBase() = default;

};//end class

}} // end namespace pressio::rom
#endif
