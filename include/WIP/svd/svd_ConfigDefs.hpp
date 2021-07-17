/*
//@HEADER
// ************************************************************************
//
// svd_ConfigDefs.hpp
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

#ifndef SVD_SVD_CONFIGDEFS_HPP_
#define SVD_SVD_CONFIGDEFS_HPP_

namespace pressio{ namespace svd{
  
// put here definitions
enum class svdType{ 
  // /*FULL SVD
  //   The matrix U'n is thus m×m, Σn is m×n diagonal, and V is n×n.
  // */
  // full, 

  // /*THIN SVD: Only the n column vectors of U corresponding to the 
  //   row vectors of V* are calculated. The remaining column 
  //   vectors of U are not calculated. This is significantly 
  //   quicker and more economical than the full SVD if n << m. 
  //   The matrix U'n is thus m×n, Σn is n×n diagonal, and V is n×n.
  // */
  // thin,

  // /*Only the r column vectors of U and r row vectors of V* 
  //   corresponding to the non-zero singular values Σr are calculated. 
  //   The remaining vectors of U and V* are not calculated. 
  //   This is quicker and more economical than the thin SVD if r << n. 
  //   The matrix Ur is thus m×r, Σr is r×r diagonal, and Vr* is r×n.*/
  // compact,


  /*Only the t column vectors of U and t row vectors of V*
    corresponding to the t largest singular values Σt are calculated. 
    The rest of the matrix is discarded. This can be much quicker 
    and more economical than the compact SVD if t≪r. The matrix Ut 
    is thus m×t, Σt is t×t diagonal, and Vt* is t×n.*/
  truncated
  
};

}} // end namespace pressio::svd
#endif  // SVD_SVD_CONFIGDEFS_HPP_
