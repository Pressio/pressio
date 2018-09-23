
#ifndef optimization_nonlinear_leastsquares_HPP_
#define optimization_nonlinear_leastsquares_HPP_

// #include "vector/core_vector_traits.hpp"
// #include "matrix/core_matrix_traits.hpp"

// #include <Eigen/Core>
// #include <unsupported/Eigen/NonLinearOptimization>
// #include <iomanip>

namespace rompp{
namespace optimization{


// template <typename matrix_type,
// 	  typename vector_type>
// typename std::enable_if< core::details::traits<matrix_type>::isEigen==1 &&
// 			 core::details::traits<vector_type>::isEigen==1
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
//   // //  jac.getNonConstRefToData() = core::details::traits<jacobian_type>::wrapped_t::Zero(jRows, jCols);
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
}//end namespace rompp
#endif

