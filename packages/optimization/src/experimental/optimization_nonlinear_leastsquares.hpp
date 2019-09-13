/*
//@HEADER
// ************************************************************************
//
// optimization_nonlinear_leastsquares.hpp
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

#ifndef optimization_nonlinear_leastsquares_HPP_
#define optimization_nonlinear_leastsquares_HPP_

// #include "vector/containers_vector_traits.hpp"
// #include "matrix/containers_matrix_traits.hpp"

// #include <Eigen/Core>
// #include <unsupported/Eigen/NonLinearOptimization>
// #include <iomanip>

namespace pressio{
namespace optimization{


// template <typename matrix_type,
// 	  typename vector_type>
// typename std::enable_if< containers::details::traits<matrix_type>::isEigen==1 &&
// 			 containers::details::traits<vector_type>::isEigen==1
// 			 >::type
// linearLstsq(const matrix_type & A,
// 	    const vector_type & b,
// 	    vector_type & x)
// {
//   // auto * eA = A.view();
//   // auto * eb = b.view();
//   // auto & ex = x.getNonConstRefToData();

//   // // ex = (*eA).bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(*eb);
//   // //ex = (*eA).colPivHouseholderQr().solve(*eb);
//   // ex = (eA->transpose() * (*eA)).ldlt().solve(eA->transpose() * *eb);
// }
// //-------------------------------------------


// template <typename vec_type>
// bool isconverged(const vec_type & y, double eps = 1e-6)
// {
//   double sum = 0.0;
//   for (decltype(y.size()) i=0; i < y.size(); i++){
//     sum += y[i]*y[i];
//   }
//   sum /= static_cast<double>(y.size());
//   sum = std::sqrt(sum);
//   return sum < eps;
// }

  
// template <typename nonlinfunctor_type,
// 	  typename state_type,
// 	  typename jacobian_type>
// void nonLinearLstsq(nonlinfunctor_type & F, state_type & y)//, int jRows, int jCols)
// {
//   // state_type res;
//   // state_type dy; dy.resize(y.size());  
//   // jacobian_type jac;
//   // //  jac.getNonConstRefToData() = containers::details::traits<jacobian_type>::wrapped_t::Zero(jRows, jCols);
//   // int maxIter = 10;
//   // // get residual and jacobian at initial guess
//   // F(y, res, jac);
  
//   // bool done = isconverged<state_type>(res, 1e-6);
//   // if (done){
//   //   return;
//   // }
//   // for (int step=0; step<maxIter; step++)
//   // {
//   //   for (decltype(y.size()) i=0; i < y.size(); i++)
//   //     dy[i] = 0.0;
//   //   linearLstsq(jac, res, dy);

//   //   for (decltype(y.size()) i=0; i < y.size(); i++){
//   //     y[i] -= dy[i];
//   //   }
//   //   F(y, res, jac);  
//   //   done = isconverged(res, 1e-6);
//   //   if (done)
//   //     	return;
//   // }//end for
  
// };


  
}//end namespace   
}//end namespace pressio
#endif

