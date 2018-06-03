

#ifndef algebra_newtonraphson_HPP
#define algebra_newtonraphson_HPP

#include "matrix/core_matrix_eigen.hpp"
#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>


namespace algebra{


template <typename matrix_type,
	  typename vector_type>
typename std::enable_if< core::details::traits<matrix_type>::isEigen==1 &&
			 core::details::traits<vector_type>::isEigen==1
			 >::type
linearSolve(const matrix_type & A,
	    const vector_type & b,
	    vector_type & x)
{
  auto * eA = A.view();
  auto * eb = b.view();
  auto & ex = x.getNonConstRefToData();

  ex = (*eA).colPivHouseholderQr().solve(*eb);
}
//-------------------------------------------


template <typename vec_type>
bool isconverged2(const vec_type & y, double eps = 1e-6)
{
  double sum = 0.0;
  for (decltype(y.size()) i=0; i < y.size(); i++){
    sum += y[i]*y[i];
  }
  sum /= static_cast<double>(y.size());
  sum = std::sqrt(sum);    
  return sum < eps;
}

  
template <typename nonlinfunctor_type,
	  typename state_type,
	  typename jacobian_type>
void newtonRaph(nonlinfunctor_type & F, state_type & y)
{
  state_type res;
  state_type dy; dy.resize(y.size());
  jacobian_type jac;
  jac.getNonConstRefToData() = core::details::traits<jacobian_type>::wrapped_t::Zero(y.size(), y.size());
  int maxIter = 100;
  // on entry, y contains the initial guess
  
  // get residual and jacobian at initial guess
  F(y, res, jac);

  // std::cout << "After EULER, Residual" << std::endl;
  // for (int i=0; i < res.size(); ++i)
  //   std::cout << res[i]  << " ";
  // std::cout << std::endl;
  // std::cout << "After EULER, jac" << std::endl;
  // std::cout << *jac.view() << std::endl;
  
  bool done = isconverged2<state_type>(res, 1e-6);
  if (done){
    return;
  }
  for (int step=0; step<maxIter; step++)
  {
    for (decltype(y.size()) i=0; i < y.size(); i++)
      dy[i] = 0.0;

    linearSolve(jac, res, dy);

    for (decltype(y.size()) i=0; i < y.size(); i++){
      y[i] -= dy[i];
    }
    F(y, res, jac);
    done = isconverged2(res, 1e-6);
    if (done)
      	return;
  }//end for

};
  
}//end namespace core
  
#endif

