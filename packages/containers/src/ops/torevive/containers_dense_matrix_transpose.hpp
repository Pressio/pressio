/*
//@HEADER
// ************************************************************************
//
// containers_dense_matrix_transpose.hpp
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

#ifndef CONTAINERS_MATRIX_OPERATIONS_DENSE_MATRIX_TRANSPOSE_HPP_
#define CONTAINERS_MATRIX_OPERATIONS_DENSE_MATRIX_TRANSPOSE_HPP_

#include "../../meta/containers_matrix_meta.hpp"
#include "../concrete/containers_matrix_dense_sharedmem_eigen.hpp"

namespace pressio{
namespace containers{
namespace mat_ops{

/*-----------------------------------------------------
  EIGEN DENSE DYNAMIC
----------------------------------------------------- */
template <typename mat_type,
	  ::pressio::mpl::enable_if_t<
	    details::traits<mat_type>::isEigen &&
	    details::traits<mat_type>::is_dense &&
	    details::traits<mat_type>::is_static==0
	    > * = nullptr>
auto transpose(const mat_type & A)
{
  mat_type res(A);
  res.data()->transposeInPlace();
  return res;
}

/*-----------------------------------------------------
  EIGEN DENSE STATIC
----------------------------------------------------- */
template <typename mat_type,
	  ::pressio::mpl::enable_if_t<
	    details::traits<mat_type>::isEigen &&
	    details::traits<mat_type>::is_dense &&
	    details::traits<mat_type>::is_static
	    > * = nullptr>
auto transpose(const mat_type & A)
{
  using scalar_t = typename details::traits<mat_type>::scalar_t;
  static constexpr int rowsIn
    = details::traits<mat_type>::wrapped_t::RowsAtCompileTime;
  static constexpr int colsIn
    = details::traits<mat_type>::wrapped_t::ColsAtCompileTime;

  using res_t = containers::Matrix<Eigen::Matrix<scalar_t, colsIn, rowsIn>>;
  res_t res;
  *res.data() = A.data()->transpose().eval();
  return res;
}
  

// /*-----------------------------------------------------
//   EPETRA DENSE MATRIX
// ----------------------------------------------------- */
// template <typename mat_type,
// 	  typename std::enable_if<
// 	    details::traits<mat_type>::isEpetra &&
// 	    details::traits<mat_type>::is_dense
// 	    >::type * = nullptr>
// auto transpose(const mat_type & A) 
// {
//   // auto & mapA = A.getDataMap();
//   auto const & mpicomm = A.commCRef();
    
//   // let's ensure that the transposed matrix cna be distributed
//   // over the current processes, so the nuber of of cols of A
//   // needs to be greater than the number of MPI processes running 
//   const int numProc = mpicomm.NumProc();
//   if( numProc > A.globalCols() ){
//     std::cout << " ERROR: trying to transpose a dense dist matrix \
// with number of columns too small to be distributed over current communicator \n";}
//   assert( numProc < A.globalRows() );

//   Epetra_Map mapT(A.globalCols(), 0, mpicomm);  
//   mat_type C( mapT, A.globalRows() );
//   C.setZero();

//   // missing impl
  
//   return C;
// }
  
  
} // end namespace mat_ops
} // end namespace containers
}//end namespace pressio
#endif
